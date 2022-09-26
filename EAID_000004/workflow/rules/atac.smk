
rule archr_build:
    """
    Build ArchR project
    """
    input:
        frags = expand("results/{sample}/fetch/fragments.tsv.gz", sample=samples_atac+samples_multiome),
        frag_ind = expand("results/{sample}/fetch/fragments.tsv.gz.tbi", sample=samples_atac+samples_multiome),
    output:
        qc_dir = directory("results_merged/atac/archr_qc"),
        project_dir = directory("results_merged/atac/archr_init"),
        arrow_dir = directory("results_merged/atac/archr_arrows")
    params:
        sample_names = samples_atac + samples_multiome,
        seed = config["archr_seed"],
        min_frags = config["archr_min_frags"],
        min_tss_enr = config["archr_min_tss_enr"]
    log:
        console = "logs/merged/atac/archr_build/console.log",
        arrow_create = "logs/merged/atac/archr_build/arrow_create.log",
        doublets = "logs/merged/atac/archr_build/doublets.log",
        save = "logs/merged/atac/archr_build/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/build_archr_project.R"

rule archr_lsi:
    """
    ArchR dimensionality reduction
    """
    input:
        project_in = "results_merged/atac/archr_init"
    output:
        project_out = directory("results_merged/atac/archr_lsi")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/merged/atac/archr_lsi/console.log",
        move = "logs/merged/atac/archr_lsi/move.log",
        lsi_atac = "logs/merged/atac/archr_lsi/lsi_atac.log",
        save = "logs/merged/atac/archr_lsi/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_lsi.R"

rule archr_harmony:
    """
    ArchR dimensionality reduction
    """
    input:
        project_in = "results_merged/atac/archr_lsi"
    output:
        project_out = directory("results_merged/atac/archr_harmony")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/merged/atac/archr_harmony/console.log",
        move = "logs/merged/atac/archr_harmony/move.log",
        umap_plot = "logs/merged/atac/archr_harmony/umap_plot.log",
        save = "logs/merged/atac/archr_harmony/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_harmony.R"

rule archr_cluster:
    """
    ArchR clustering
    """
    input:
        project_in = "results_merged/atac/archr_harmony"
    output:
        project_out = directory("results_merged/atac/archr_clustered")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/merged/atac/archr_cluster/console.log",
        move = "logs/merged/atac/archr_cluster/move.log",
        cluster_atac = "logs/merged/atac/archr_cluster/cluster_atac.log",
        umap_plot = "logs/merged/atac/archr_cluster/umap_plot.log",
        save = "logs/merged/atac/archr_cluster/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_cluster.R"

rule archr_linkage:
    """
    ArchR unconstrained cluster linkage using reference
    """
    input:
        project_in = "results_merged/atac/archr_clustered",
        seurat_data = "reference/seurat_build_reference_log1p/proj.rds"
    output:
        project_out = directory("results_merged/atac/archr_linkage"),
        labels = "results_merged/atac/archr_label_data.tsv"
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/merged/atac/archr_linkage/console.log",
        move = "logs/merged/atac/archr_linkage/move.log",
        umap_plot = "logs/merged/atac/archr_linkage/umap_plot.log",
        linkage = "logs/merged/atac/archr_linkage/linkage.log",
        save = "logs/merged/atac/archr_linkage/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_linkage.R"

rule archr_write_qc:
    """
    ArchR write sample QC data
    """
    input:
        qc_dir = "results_merged/atac/archr_qc"
    output:
        out_dir = directory("results_merged/atac/archr_qc_parsed")
    log:
        console = "logs/merged/atac/archr_write_qc/console.log"
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_write_qc.R"

rule archr_write_data:
    """
    ArchR write cell data
    """
    input:
        project_in = "results_merged/atac/archr_linkage"
    output:
        markers = directory("results_merged/atac/archr_write_data/markers"),
        emb_coords = "results_merged/atac/archr_write_data/emb_coords.tsv",
        umap_coords = "results_merged/atac/archr_write_data/umap_coords.tsv",
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/merged/atac/archr_write_data/console.log",
        markers = "logs/merged/atac/archr_write_data/markers.log",
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_write_data.R"

rule export_atac_embeddings:
    """
    Export ATAC embeddings
    """
    input:
        emb = "results_merged/atac/archr_write_data/emb_coords.tsv",
        umap = "results_merged/atac/archr_write_data/umap_coords.tsv"
    output:
        emb = "export/atac/embeddings/harmony.tsv.gz",
        umap = "export/atac/embeddings/umap.tsv.gz"
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_atac_embeddings.py"

rule export_atac_labels:
    """
    Export ATAC cell types
    """
    input:
        "results_merged/atac/archr_label_data.tsv"
    output:
        "export/atac/labels/cell_types.tsv.gz"
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_atac_labels.py"

rule export_atac_markers:
    """
    Export ATAC markers
    """
    input:
        markers = "results_merged/atac/archr_write_data/markers",
        gtf = "results_merged/fetch/GRCh38.gtf.gz"
    output:
        directory("export/atac/markers")
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_atac_markers.py"

rule export_atac_metadata:
    """
    Export ATAC metadata
    """
    input:
        metadata = "results_merged/atac/archr_qc_parsed",
        final_data = "results_merged/atac/archr_label_data.tsv"
    output:
        "export/atac/metadata.tsv.gz",
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_atac_metadata.py"

rule export_atac_figures:
    """
    Export ATAC figures
    """
    input:
        "results_merged/atac/archr_linkage"
    output:
        scratch = directory("results_merged/atac/export_figures"),
        tarball = "export/atac/figures.tar.gz"
    params:
        readme = workflow.source_path("../resources/atac_figures_readme.txt")
    conda:
        "../envs/fetch.yaml"
    shell:
        "mkdir -p {output.scratch}; "
        "cp {input}/Plots/umap_full_label.pdf {output.scratch}/umap_labels.pdf; "
        "cp {input}/Plots/umap_harmony_datasets.pdf {output.scratch}/umap_samples.pdf; "
        "cp {params.readme} {output.scratch}/README.txt; "
        "tar -zcvf {output.tarball} {output.scratch}"

rule export_atac_dataset_names:
    """
    Export ATAC dataset names used in analysis
    """
    output:
        "export/atac/datasets.txt"
    params:
        datasets = [i.split("_")[-1].split("-")[0] for i in samples_atac+samples_multiome]
    conda:
        "../envs/fetch.yaml"
    shell:
        "echo {params.datasets} | tr ' ' '\\n' > {output}"