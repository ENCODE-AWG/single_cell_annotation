
scale_fill_rainbow <- function(saturation=4, ...) {
  rainbow_palette <- RColorBrewer::brewer.pal(10, "Spectral") %>% multiomeUtils:::saturate(saturation) %>% rev()
  scale_fill_gradientn(colors=rainbow_palette, ...)
}

cutoff_plot_hexbin <- function(data, x, y, x_cutoff, y_cutoff) {
  #browser()
  quadrant_counts <- data %>% 
    summarize(
      top_right = sum({{x}} > x_cutoff & {{y}} > y_cutoff),
      top_left = sum({{x}} <= x_cutoff & {{y}} > y_cutoff),
      bot_right = sum({{x}} > x_cutoff & {{y}} <= y_cutoff),
      bot_left = sum({{x}} <= x_cutoff & {{y}} <= y_cutoff)
    ) %>%
    mutate(across(everything(), scales::comma))
  
  quadrant_label <- tribble(
    ~x, ~y, ~hjust, ~vjust, ~text,
    Inf,  Inf,  1.5,  2, quadrant_counts$top_right,
    -Inf,  Inf,  -0.5,  2, quadrant_counts$top_left,
    Inf, -Inf,  1.5, -2, quadrant_counts$bot_right,
    -Inf, -Inf,  -0.5, -2, quadrant_counts$bot_left
  )
  
  ggplot(data, aes({{x}}, {{y}})) +
    stat_bin_hex(bins = 100) +
    scale_fill_rainbow(trans="log10") +
    geom_hline(yintercept=y_cutoff, linetype="dashed") +
    geom_vline(xintercept=x_cutoff, linetype="dashed") +
    geom_text(data=quadrant_label, mapping=aes(x=x, y=y, hjust=hjust, vjust=vjust, label=text)) +
    guides(fill="none")
}

#' Make a knee plot for displaying results from read_count_cutoff
#' @param read_counts Vector of read counts per cell
#' @param cutoff Read cutoff calculated from read_count_cutoff
#' @export
read_count_cutoff_plot <- function(read_counts, cutoff) {
  # Downsample the total cells considered while fitting spline.
  # Keep ~1k cells per order of magnitude
  ranks <- unique(floor(10^(1:(1000*ceiling(log10(length(read_counts))))/1000)))
  ranks <- ranks[ranks < length(read_counts)]
  
  reads <- sort(read_counts, decreasing=TRUE)[ranks]
  min_rank <- sum(read_counts > cutoff)
  
  ggplot(NULL, aes(log10(ranks), reads)) +
    geom_point() +
    scale_x_continuous(labels=scales::label_math()) +
    scale_y_log10() +
    geom_hline(yintercept=cutoff, linetype="dashed") +
    labs(x = "Barcode Rank", y="Reads",
         subtitle=sprintf("%s cells with > %s reads", 
                          scales::comma(min_rank), 
                          scales::comma(floor(cutoff))))
}

#' Load a gtf file into a GRanges object using data.table::fread.
#' I'm sure there's other ways of doing this, but this is simple and works for me.
#' @param path File path of gtf
#' @param attributes Attribute names to be parsed out into separate metadata columns.
#'    These show up in the last column of the gtf as attribute_name "value";, e.g. transcript_id "ENST00000619216.1";
#' @param keep_attribute_column Boolean for whether to preserve the text list of attributes as a metadata column.
#'    (Useful if additional parsing of attributes is needed)
#' @param ... Additional arguments passed to fread. Of particular use is skip, for skipping commented lines at the start.
#'    Set it to the text at the start of the first line to read. skip="chr" is often appropriate
#' @return GRanges object with one entry per line of the input GTF.
read_gtf <- function(path, attributes, keep_attribute_colum=FALSE, ...) {
  gtf_colnames <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
  
  ret <- data.table::fread(path, col.names=gtf_colnames, ...)
  for (a in attributes) {
    ret[[a]] <- str_match(ret[["attributes"]], sprintf('%s "([^"]*)"', a))[,2]
  }
  ret <- GenomicRanges::makeGRangesFromDataFrame(ret, keep.extra.columns = TRUE)
  if(!keep_attribute_colum)
    ret[["attributes"]] <- NULL
  return(ret)
}



#' Get a knn matrix given SVD matrix for reference and query
#' 
#' @param data cell x dims matrix for reference dataset
#' @param query cell x dims matrix for query dataset (optional)
#' @param k number of neighbors to calculate
#' @param threads Number of threads to use. Note that result is non-deterministic
#'          if threads > 1
#' @param ef ef parameter for RccppHNSW::hnsw_search. Increase for slower search but
#'          improved accuracy
#' @param verbose whether to print progress information during search
#' @return List of 2 matrices -- idx for cell x K neighbor indices,
#'         dist for cell x K neighbor distances.
#'         If no query is given, nearest neighbors are found mapping
#'         the data matrix to itself, prohibiting self-neighbors
#' @export
get_knn <- function(data, query=NULL, k = 10, metric="euclidean", verbose=TRUE, threads=1, ef=100) {
  index <- RcppHNSW::hnsw_build(
    data, 
    distance=metric, 
    verbose=verbose,
    n_threads=threads
  )
  if (is.null(query)) {
      res <- RcppHNSW::hnsw_search(
        data,
        index,
        k+1,
        ef=ef,
        verbose=verbose,
        n_threads = threads
      )
      missed_searches <- sum(res$idx[,1] != seq_len(nrow(data)))
      if (any(res$idx[,1] != seq_len(nrow(data)))) {
        stop(sprintf("KNN search didn't find self-neighbor for %d datapoints", missed_searches))
      }
      return(list(
        idx = res$idx[,-1,drop=FALSE],
        dist = res$dist[,-1,drop=FALSE]
      ))
  } else {
      res <- RcppHNSW::hnsw_search(
        query,
        index,
        k,
        ef=ef,
        verbose=verbose,
        n_threads = threads
      )
      return(res)
  }
}


#' Internal worker function for label transfer and quality estimation
#' @param knn KNN object as returned from `get_knn`.
#' @param labels Vector of label values to transfer. Should be type integer with minimum >= 1
#'               (convert from strings or factors before calling this function)
#' @return Matrix of cells x labels
get_label_scores <- function(knn, labels) {
  stopifnot(is.integer(labels))
  stopifnot(min(labels) >= 1)

  nearest_labels <- matrix(labels[knn$idx], nrow=nrow(knn$idx))
  
  weights <- 1 - knn$dist / knn$dist[,ncol(knn$dist)]
  weights <- weights / rowSums(weights)
 
  label_scores <- map(
    seq_len(max(labels)),
    ~ rowSums(weights * (nearest_labels == .x))                    
  ) %>% do.call(cbind, .)
  
  return(label_scores)
}

#' Transfer labels using the most common value from the K nearest neighbors
#' 
#' Neighbors are weighted similarly to Stuart et al. (2019), by 1 - (dist_i/dist_k)
#' @param knn KNN object as returned from `get_knn`.
#' @param labels Vector of label values to transfer
#' @export
transfer_labels <- function(knn, labels) {
  label_values <- unique(labels)
  label_indices <- match(labels, label_values)
  
  label_scores <- get_label_scores(knn, label_indices)
  best_label <- max.col(label_scores, ties.method = "first")
  
  label_values[best_label]
}