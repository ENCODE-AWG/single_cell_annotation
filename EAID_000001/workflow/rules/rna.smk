rule seurat_build_reference:
    """
    Build Seurat reference dataset
    """
    input:
        mat = "reference/fetch/matrix.mtx",
        features = "reference/fetch/features.tsv",
        cells = "reference/fetch/barcodes.tsv",
        metadata = "reference/fetch/metadata.csv"
    output:
        project_out = "reference/seurat_build_reference/proj.rds",
        qc_violin = "reference/seurat_build_reference/qc_violin.pdf",
        qc_scatter = "reference/seurat_build_reference/qc_scatter.pdf",
        umap =  "reference/seurat_build_reference/umap.pdf"
    params:
        seed = config["seurat_seed"]
    log:
        console = "logs/reference/seurat_build_reference/seurat_build_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/build_seurat_reference.R"

rule seurat_build_reference_log1p:
    """
    Build Seurat reference dataset (no SCT)
    """
    input:
        mat = "reference/fetch/matrix.mtx",
        features = "reference/fetch/features.tsv",
        cells = "reference/fetch/barcodes.tsv",
        metadata = "reference/fetch/metadata.csv"
    output:
        project_out = "reference/seurat_build_reference_log1p/proj.rds",
        qc_violin = "reference/seurat_build_reference_log1p/qc_violin.pdf",
        qc_scatter = "reference/seurat_build_reference_log1p/qc_scatter.pdf",
        umap =  "reference/seurat_build_reference_log1p/umap.pdf"
    params:
        seed = config["seurat_seed"]
    log:
        console = "logs/reference/seurat_build_reference_log1p/seurat_build_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/build_seurat_reference_log1p.R"

rule seurat_build_rna:
    """
    Build Seurat project
    """
    input:
        mat = "results/{sample}/fetch/matrix.mtx",
        features = "results/{sample}/fetch/features.tsv",
        cells = "results/{sample}/fetch/barcodes.tsv"
    output:
        project_out = "results/{sample}/rna/seurat_build_rna/proj.rds",
        metadata = "results/{sample}/rna/seurat_build_rna/metadata.tsv",
        qc_violin = "results/{sample}/rna/seurat_build_rna/qc_violin.pdf",
        qc_scatter = "results/{sample}/rna/seurat_build_rna/qc_scatter.pdf",
    params:
        sample_name = lambda w: w.sample,
        seed = config["seurat_seed"],
        min_count_rna = config["seurat_min_count"],
        max_pct_mito_rna = config["seurat_max_pct_mito"]
    log:
        console = "logs/{sample}/rna/seurat_build_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/build_seurat_rna.R"

rule seurat_doublets_rna:
    """
    Filter doublets
    """
    input:
        project_in = "results/{sample}/rna/seurat_build_rna/proj.rds"
    output:
        project_out_all = "results/{sample}/rna/seurat_doublets_rna/proj_all.rds",
        project_out_filtered = "results/{sample}/rna/seurat_doublets_rna/proj_filtered.rds"
    params:
        seed = config["seurat_seed"],
        doublet_rate = config["doublet_formation_rate"]
    log:
        console = "logs/{sample}/rna/seurat_doublets_rna/console.log"
    threads:
        config["max_threads_per_rule"]
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_doublets_rna.R"

rule seurat_label_rna:
    """
    Seurat RNA reference labeling
    """
    input:
        project_rna = "results/{sample}/rna/seurat_doublets_rna/proj_filtered.rds",
        project_ref = "reference/seurat_build_reference/proj.rds"
    output:
        project_out = "results/{sample}/rna/seurat_label_rna/proj.rds",
        umap = "results/{sample}/rna/seurat_label_rna/umap.pdf"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/{sample}/rna/seurat_label_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_label_rna.R"

rule seurat_merge_rna:
    """
    Merge RNA samples
    """
    input:
        projects_in = expand("results/{sample}/rna/seurat_label_rna/proj.rds", sample=samples_rna+samples_multiome)
    output:
        project_out = "results_merged/rna/seurat_merge_rna/proj.rds",
        umap_dataset_pre_harmony = "results_merged/rna/seurat_merge_rna/umap_dataset_pre_harmony.pdf",
        umap_pre_harmony = "results_merged/rna/seurat_merge_rna/umap_pre_harmony.pdf",
        umap_dataset = "results_merged/rna/seurat_merge_rna/umap_dataset.pdf",
        umap = "results_merged/rna/seurat_merge_rna/umap.pdf",
    params:
        seed = config["seurat_seed"],
        samples = samples_rna + samples_multiome,
    log:
        console = "logs/merged/rna/seurat_merge_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_merge_rna.R"

rule seurat_cluster_rna:
    """
    Seurat RNA clustering
    """
    input:
        project_in = "results_merged/rna/seurat_merge_rna/proj.rds",
        # project_ref = "reference/seurat_build_reference/proj.rds"
    output:
        project_out = "results_merged/rna/seurat_cluster_rna/proj.rds",
        umap = "results_merged/rna/seurat_cluster_rna/umap_clusters.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/rna/seurat_cluster_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_cluster_rna.R"

rule seurat_name_rna:
    """
    Seurat RNA cluster naming
    """
    input:
        project_in = "results_merged/rna/seurat_cluster_rna/proj.rds",
    output:
        project_out = "results_merged/rna/seurat_name_rna/proj.rds",
        umap = "results_merged/rna/seurat_name_rna/umap_clusters.pdf",
        umap_qc = "results_merged/rna/seurat_name_rna/umap_qc.pdf",
        mat = "results_merged/rna/seurat_name_rna/confusion_mat.pdf",
        metadata = "results_merged/rna/seurat_name_rna/metadata.tsv",
    params:
        seed = config["seurat_seed"],
        # cluster_names = rna_cluster_names
    log:
        console = "logs/merged/rna/seurat_name_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_name_rna.R"

rule seurat_write_embeddings:
    """
    Seurat save RNA embeddings
    """
    input:
        project_in = "results_merged/rna/seurat_name_rna/proj.rds",
    output:
        emb_coords = "results_merged/rna/seurat_write_embeddings/emb_coords.tsv",
        umap_coords = "results_merged/rna/seurat_write_embeddings/umap_coords.tsv",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/rna/seurat_write_embeddings/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_write_embeddings.R"

rule seurat_write_markers:
    """
    Seurat save RNA marker genes
    """
    input:
        project_in = "results_merged/rna/seurat_name_rna/proj.rds",
    output:
        markers = directory("results_merged/rna/seurat_write_markers/markers"),
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/rna/seurat_write_markers/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_write_markers.R"

rule export_rna_embeddings:
    """
    Export RNA embeddings
    """
    input:
        emb = "results_merged/rna/seurat_write_embeddings/emb_coords.tsv",
        umap = "results_merged/rna/seurat_write_embeddings/umap_coords.tsv",
    output:
        emb = "export/rna/embeddings/harmony.tsv.gz",
        umap = "export/rna/embeddings/umap.tsv.gz",
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_rna_embeddings.py"

rule export_rna_labels:
    """
    Export RNA cell types
    """
    input:
        "results_merged/rna/seurat_name_rna/metadata.tsv"
    output:
        "export/rna/labels/cell_types.tsv.gz",
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_rna_labels.py"

rule export_rna_markers:
    """
    Export RNA markers
    """
    input:
        markers = "results_merged/rna/seurat_write_markers/markers",
        genes = expand("results/{sample}/fetch/features.tsv", sample=samples_rna+samples_multiome)
    output:
        directory("export/rna/markers")
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_rna_markers.py"

rule export_rna_metadata:
    """
    Export RNA metadata
    """
    input:
        metadata = expand("results/{sample}/rna/seurat_build_rna/metadata.tsv", sample=samples_rna+samples_multiome),
        final_data = "results_merged/rna/seurat_name_rna/metadata.tsv"
    output:
        "export/rna/metadata.tsv.gz",
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_rna_metadata.py"

rule export_rna_figures:
    """
    Export RNA figures
    """
    input:
        umap_labels = "results_merged/rna/seurat_name_rna/umap_clusters.pdf",
        umap_samples = "results_merged/rna/seurat_merge_rna/umap_dataset.pdf"
    output:
        scratch = directory("results_merged/rna/export_figures"),
        tarball = "export/rna/figures.tar.gz"
    params:
        readme = workflow.source_path("../resources/rna_figures_readme.txt")
    conda:
        "../envs/fetch.yaml"
    shell:
        "mkdir -p {output.scratch}; "
        "cp {input.umap_labels} {output.scratch}/umap_labels.pdf; "
        "cp {input.umap_samples} {output.scratch}/umap_samples.pdf; "
        "cp {params.readme} {output.scratch}/README.txt; "
        "tar -zcvf {output.tarball} {output.scratch}"

rule export_rna_dataset_names:
    """
    Export RNA dataset names used in analysis
    """
    output:
        "export/rna/datasets.txt"
    params:
        datasets = [i.split("_")[0].split("-")[0] for i in samples_rna+samples_multiome]
    conda:
        "../envs/fetch.yaml"
    shell:
        "echo {params.datasets} | tr ' ' '\\n' > {output}"