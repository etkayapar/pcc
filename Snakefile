include: "common_utils.smk"

rule all:
    input:
        "outlier_detection/outlier_genes.txt"

rule align_aa:
    input:
        "genewise_fastas/{gene}.faa"
    output:
        "mafft_output/{gene}_aligned.faa"
    threads: 8
    conda:
        "envs/mafft.yaml"
    shell:
        "linsi --thread {threads} {input} > {output}"

rule backtranslate:
    input:
        nt="genewise_fastas/{gene}.fa",
        aa_msa="mafft_output/{gene}_aligned.faa",
        pal2nal_path=workflow.basedir+"/utils/extract-buscos/pal2nal.py"
    output:
        "mafft_output/{gene}_aligned.fa"
    log:
        workflow.basedir+"/logs/backtranslate/{gene}_backtranslate.log"
    conda:
        "envs/phylo_scripts_python.yaml"
    shell:
        """
        ln -r -s {input.nt} mafft_output/{wildcards.gene}.fa
        cd mafft_output
        (python3 {input.pal2nal_path} {wildcards.gene}) 2>tomim.err
        """

rule clean_all_gap_seqs:
    input:
        rules.backtranslate.output
    output:
        "before_trimal/{gene}.fa"
    shell:
        "utils/phylo_scripts/cleanAllGaps {input} > {output}"


rule infer_gene_trees_before_trimal:
    input:
        "before_trimal/{gene}.fa"
    output:
        treefile="gene_trees/before_trimal/{gene}/{gene}.treefile",
        treedir=directory("gene_trees/before_trimal/{gene}")
    threads: 4
    conda:
        "envs/iqtree.yaml"
    shell:
        """
        iqtree2 -s {input} \
                -m MFP -mset GTR -mrate I+R \
                -T {threads} --prefix {wildcards.gene} -st DNA --keep-ident
        mv {wildcards.gene}.* {output.treedir}/
        """
# mkdir -p {wildcards.gene}
# mv {wildcards.gene}.* {wildcards.gene}/
# mv {wildcards.gene} {output.treedir}

rule collect_gene_trees:
    input:
        expand("gene_trees/before_trimal/{gene}/{gene}.treefile", gene=genes)
    output:
        trees="outlier_detection/all_genes.treefile",
        gene_names="outlier_detection/all_genes_names.txt"
    shell:
        """
        ls {input} | cut -d '/' -f3 > {output.gene_names}
        cat {input} > {output.trees}
        """

rule detect_outliers:
    input:
        treefile=rules.collect_gene_trees.output.trees,
        gene_names=rules.collect_gene_trees.output.gene_names
    output:
        saved_genes_plot="outlier_detection/saved_genes.pdf",
        outlier_genes_plot="outlier_detection/outlier_genes.pdf",
        outlier_genes_list="outlier_detection/outlier_genes.txt",
        kept_taxa_dir=directory("outlier_detection/saved_genes_kept_taxa")
    conda:
        "envs/detect_outliers.yaml"
    params:
        long_branch_threshold=0.1,
        taxa_threshold=0.7
    shell:
        """
        Rscript utils/gene_trees.R {params.long_branch_threshold} {params.taxa_threshold} \
        {input.treefile} {input.gene_names} outlier_detection 
        """
