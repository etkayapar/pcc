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
gene_names_path=snakemake.input[1]
outdir=os.path.dirname(snakemake.output[0])

## Run TreeShrink
#run_treeshrink.py -t {input.treefile} -m "per-gene" -o {output.removed_taxa_dir}
sp.run(['run_treeshrink.py', '-t', snakemake.input[0], '-m', treeshrink_mode,
        '-o', outdir])
    
