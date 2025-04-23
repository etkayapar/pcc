"""
Python version of the detect_outliers_before_trimal and
detect_outliers_after_trimal rules
"""

import os
import subprocess as sp

taxon_threshold=snakemake.params[0]
pipeline_stage=snakemake.params[1]
treeshrink_mode=snakemake.params[2]
long_branch_threshold=snakemake.params[3]
r_script_path="utils/detect_paralogs.R"
treeshrink_output_trees_path=snakemake.input[2]
gene_names_path=snakemake.input[1]
outdir=os.path.dirname(treeshrink_output_trees_path)

## Touch
open(snakemake.output[0],'a').close()

## Calc ngenes
with open(snakemake.input[1], 'r') as f:
    genenames = [x.strip() for x in f.readlines()]
removed_taxa_path = os.path.join(outdir, "output.txt")
with open(removed_taxa_path, 'r') as f:
    removed_taxa = [x.strip().split() for x in f.readlines()]

## Here will come running and parsing the gene_trees Rscript
sp.run(['Rscript', r_script_path, str(long_branch_threshold),
        "0.55", treeshrink_output_trees_path,
        gene_names_path, outdir])

for gene,genename in enumerate(genenames):
    outlier_genes_path = snakemake.output[0]
    outlier_genes_file = open(outlier_genes_path, "w")
    additional_removed_taxa_path = os.path.join(outdir, genename+"_removed_taxa_r.txt")
    if os.path.exists(additional_removed_taxa_path):
        with open(additional_removed_taxa_path, "r") as f:
            additional_taxa = [x.strip() for x in f.readlines()]
        removed_taxa[gene].extend(additional_taxa)
    this_gene_removed_taxa = removed_taxa[gene]
    with open(f"output/{pipeline_stage}/gene_tree_input/{genename}.fa") as f:
        fasta = f.readlines()
    original_ntax = len([x for x in fasta if x.startswith(">")])
    new_ntax = original_ntax - len(this_gene_removed_taxa)
    taxon_retaining_pct = new_ntax / original_ntax * 100
    print(f"{genename}\t{original_ntax}\t{new_ntax}\t{round(taxon_retaining_pct,2)}")
    if taxon_retaining_pct == 100:
        continue
    if taxon_retaining_pct != 0 and taxon_retaining_pct >= taxon_threshold:
        with open(os.path.join(outdir, f"{genename}_removed_taxa.txt"), 'w') as f:
            for taxon in this_gene_removed_taxa:
                f.write(taxon+'\n')
    else:
        outlier_genes_file.write(genename+'\n')

outlier_genes_file.close()
                
    
        
