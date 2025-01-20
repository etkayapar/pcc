"""
Python version of the detect_outliers_before_trimal and
detect_outliers_after_trimal rules
"""

import os
import subprocess as sp

pipeline_stage=snakemake.params[1]

## Touch
open(snakemake.output[0],'a').close()

## Run TreeShrink
#run_treeshrink.py -t {input.treefile} -m "per-gene" -o {output.removed_taxa_dir}
sp.run(['run_treeshrink.py', '-t', snakemake.input[0], '-m', 'auto', '-o', snakemake.output[1]])
## Calc ngenes
with open(snakemake.input[1], 'r') as f:
    genenames = [x.strip() for x in f.readlines()]
removed_taxa_path = os.path.join(snakemake.output[1], "output.txt")
with open(removed_taxa_path, 'r') as f:
    removed_taxa = [x.strip().split() for x in f.readlines()]

for gene,genename in enumerate(genenames):
    outlier_genes_path = snakemake.output[0]
    outlier_genes_file = open(outlier_genes_path, "w")
    # breakpoint()
    this_gene_removed_taxa = removed_taxa[gene]
    with open(f"output/{pipeline_stage}/gene_tree_input/{genename}.fa") as f:
        fasta = f.readlines()
    original_ntax = len([x for x in fasta if x.startswith(">")])
    new_ntax = original_ntax - len(this_gene_removed_taxa)
    taxon_retaining_pct = new_ntax / original_ntax * 100
    # breakpoint()
    if taxon_retaining_pct == 100:
        continue
    if taxon_retaining_pct != 0 and taxon_retaining_pct >= snakemake.params[0]:
        with open(os.path.join(snakemake.output[1], f"{genename}_removed_taxa.txt"), 'w') as f:
            for taxon in this_gene_removed_taxa:
                f.write(taxon+'\n')
    else:
        outlier_genes_file.write(genename+'\n')
                
    
        
