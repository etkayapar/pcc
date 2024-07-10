rule init_after_trimal:
    input:
        "before_trimal/outlier_detection/realignment/{gene}_aligned.fa"
    output:
        "after_trimal/gene_tree_input/{gene}.fa"
    shell:
        """
        utils/phylo_scripts/cleanAllGaps {input} > {output}
        """

rule detect_outliers_after_trimal:
    input:
        treefile="after_trimal/outlier_detection/all_genes.treefile",
        gene_names="after_trimal/outlier_detection/all_genes_names.txt"
    output:
        saved_genes_plot  ="after_trimal/outlier_detection/saved_genes.pdf",
        outlier_genes_plot="after_trimal/outlier_detection/outlier_genes.pdf",
        outlier_genes_list="after_trimal/outlier_detection/outlier_genes.txt",
        kept_taxa_dir=directory("after_trimal/outlier_detection/saved_genes_kept_taxa")
    conda:
        "../envs/detect_outliers.yaml"
    params:
        long_branch_threshold=0.1,
        taxa_threshold=0.7
    shell:
        f"""
        Rscript utils/gene_trees.R {{params.long_branch_threshold}} {{params.taxa_threshold}} {{input.treefile}} {{input.gene_names}} after_trimal/outlier_detection > gt.log 2>gt.err
        """

checkpoint process_outliers_after_trimal:
    input:
        aln_dir="after_trimal/gene_tree_input",
        keep_taxa_path="utils/phylo_scripts/keep_taxa.awk"
    output:
        genelist="after_trimal/outlier_detection/final_output/genelist.txt",
        d=directory("after_trimal/outlier_detection/final_output")
        #"{stage}/outlier_detection/final_output/{gene}.fa"
    params:
        #kept_taxa_path="after_trimal/outlier_detection/saved_genes_kept_taxa/{{gene}}_kept_taxa.txt",
        outlier_genes_path="after_trimal/outlier_detection/outlier_genes.txt"
    log:
        workflow.basedir+"/logs/after_trimal/process_outliers.log"
    shell:
        """
        set +o pipefail
        (for gene in `ls {input.aln_dir}/ | grep -E "*.fa$" | cut -d. -f 1`
        do
            kept_taxa_path="after_trimal/outlier_detection/saved_genes_kept_taxa/${{gene}}_kept.taxa.txt"
            if [[ -f $kept_taxa_path ]]
            then
            {input.keep_taxa_path} -v taxafile=${{kept_taxa_path}} {input.aln_dir}/${{gene}}.fa > {output.d}/${{gene}}.fa
            echo ${{gene}} >> {output.genelist}
            elif grep -qw ${{gene}} {params.outlier_genes_path}
            then
            true
            else
            ln -sr {input.aln_dir}/${{gene}}.fa {output.d}
            echo ${{gene}} >> {output.genelist}
            fi
        done) 2>{log}
        """

rule backtranslate_final:
    input:
        nt="after_trimal/outlier_detection/realignment/{gene}.fa",
        aa_msa="after_trimal/outlier_detection/realignment/{gene}_aligned.faa",
        pal2nal_path=workflow.basedir+"/utils/extract-buscos/pal2nal.py"
    output:
        nt_aln="after_trimal/outlier_detection/realignment/{gene}_aligned.fa"
    log:
        workflow.basedir+"/logs/after_trimal/backtranslate_final/{gene}_backtranslate.log"
    conda:
        "../envs/phylo_scripts_python.yaml"
    shell:
        """
        cd $(dirname {output.nt_aln})
        (python3 {input.pal2nal_path} {wildcards.gene}) 2> {log}
        """

rule init_concatenate:
    input:
        get_gene_list_to_concatenate
    output:
        d=directory("genes_to_concat/"),
        okf="genes_to_concat/OK"
    shell:
        """
        mkdir -p {output.d}
        cp {input} {output.d}
        touch {output.okf}
        """


rule concatenate:
    input:
        okf="genes_to_concat/OK"
    output:
        data="supermatrix.phy",
        part="supermatrix.nex"
    conda:
        "../envs/phylo_scripts_python.yaml"
    shell:
        """
        rm {input.okf}
        concat-aln $(dirname {input.okf}) supermatrix DNA
        """
