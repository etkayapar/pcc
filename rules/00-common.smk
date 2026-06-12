import os
import pandas as pd

configfile: "config.yaml"

gene_table = pd.read_table(config["gene_table"], index_col=False, dtype=str)
genes = gene_table["gene"]

def get_gene_list_to_infer_tree_after(wildcards):
    gene_tree_dir = checkpoints.process_outliers_before_trimal.get(**wildcards).output[1]
    genes_list = expand("output/after_trimal/gene_trees/{gene}/{gene}.treefile",gene=glob_wildcards(os.path.join(gene_tree_dir, '{gene}.fa')).gene)

    return genes_list

def get_filtered_genes_after_first_pass(wildcards):
    with checkpoints.process_outliers_before_trimal.get(**wildcards).output[0].open() as f:
        genes = [line.strip() for line in f]
    return genes

def get_gene_list_to_concatenate(wildcards):
    gene_tree_dir = checkpoints.process_outliers_after_trimal.get(**wildcards).output[1]
    genes_list = expand("output/after_trimal/outlier_detection/realignment/{gene}_aligned.fa",gene=glob_wildcards(os.path.join(gene_tree_dir, '{gene}.fa')).gene)

    return genes_list

def get_filtered_genes_final(wildcards):
    with checkpoints.process_outliers_after_trimal.get(**wildcards).output[0].open() as f:
        genes = [line.strip() for line in f]
    return genes

def get_aln_params(wildcards, input):
    aligner = config["params"]["align_aa"]["aligner"]
    if aligner != "auto":
        if aligner not in ["einsi", "linsi", "ginsi", "fftns", "fftnsi"]:
            raise ValueError(
                "Configured aligner not supported. Choose from 'auto', 'einsi',"
                "'linsi', 'ginsi', 'fftns', or 'fftnsi'."
            )
        return aligner
    auto_criterion = config["params"]["align_aa"]["auto_criterion"]
    if auto_criterion != "filesize":
        if auto_criterion != "mafft":
            raise ValueError("The 'auto_criterion' should be either 'mafft' or 'filesize'")
        return "mafft --auto"
    large_gene_threshold_mb = config["params"]["align_aa"]["large_gene_threshold_mb"]
    small_gene_aligner = config["params"]["align_aa"]["small_gene_aligner"]
    large_gene_aligner = config["params"]["align_aa"]["large_gene_aligner"]
    file_size_mb = os.path.getsize(input[0]) / (1024 * 1024)
    return large_gene_aligner if file_size_mb > large_gene_threshold_mb else small_gene_aligner
