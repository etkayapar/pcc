import os
import pandas as pd

configfile: "config.yaml"

gene_table = pd.read_table(config["gene_table"], index_col=False, dtype=str)
genes = gene_table["gene"]

def get_gene_list_to_infer_tree_after(wildcards):
    gene_tree_dir = checkpoints.process_outliers_before_trimal.get(**wildcards).output[1]
    genes_list = expand("output/after_trimal/gene_trees/{gene}/{gene}.treefile",gene=glob_wildcards(os.path.join(gene_tree_dir, '{gene}.fa')).gene)

    return genes_list

def get_gene_list_to_concatenate(wildcards):
    gene_tree_dir = checkpoints.process_outliers_after_trimal.get(**wildcards).output[1]
    genes_list = expand("output/after_trimal/outlier_detection/realignment/{gene}_aligned.fa",gene=glob_wildcards(os.path.join(gene_tree_dir, '{gene}.fa')).gene)

    return genes_list

