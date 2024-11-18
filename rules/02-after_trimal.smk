rule init_after_trimal:
    input:
        "output/before_trimal/outlier_detection/realignment/{gene}_aligned.fa"
    output:
        "output/after_trimal/gene_tree_input/{gene}.fa"
    shell:
        """
        utils/phylo_scripts/cleanAllGaps {input} > {output}
        """

rule detect_outliers_after_trimal:
    input:
        treefile="output/after_trimal/outlier_detection/all_genes.treefile",
        gene_names="output/after_trimal/outlier_detection/all_genes_names.txt"
    output:
        saved_genes_plot  ="output/after_trimal/outlier_detection/saved_genes.pdf",
        outlier_genes_plot="output/after_trimal/outlier_detection/outlier_genes.pdf",
        outlier_genes_list="output/after_trimal/outlier_detection/outlier_genes.txt",
        kept_taxa_dir=directory("output/after_trimal/outlier_detection/saved_genes_kept_taxa")
    conda:
        "../envs/detect_outliers.yaml"
    params:
        long_branch_threshold=config["params"]["detect_outliers"]["long_branch_threshold"],
        taxa_threshold=config["params"]["detect_outliers"]["taxa_threshold"]
    shell:
        f"""
        Rscript utils/gene_trees.R {{params.long_branch_threshold}} {{params.taxa_threshold}} {{input.treefile}} {{input.gene_names}} output/after_trimal/outlier_detection > gt.log 2>gt.err
        """

checkpoint process_outliers_after_trimal:
    input:
        aln_dir="output/after_trimal/gene_tree_input",
        keep_taxa_path="utils/phylo_scripts/keep_taxa.awk"
    output:
        genelist="output/after_trimal/outlier_detection/final_output/genelist.txt",
        d=directory("output/after_trimal/outlier_detection/final_output")
        #"{stage}/outlier_detection/final_output/{gene}.fa"
    params:
        #kept_taxa_path="output/after_trimal/outlier_detection/saved_genes_kept_taxa/{{gene}}_kept_taxa.txt",
        outlier_genes_path="output/after_trimal/outlier_detection/outlier_genes.txt"
    log:
        workflow.basedir+"/logs/after_trimal/process_outliers.log"
    shell:
        """
        set +o pipefail
        (for gene in `ls {input.aln_dir}/ | grep -E "*.fa$" | cut -d. -f 1`
        do
            kept_taxa_path="output/after_trimal/outlier_detection/saved_genes_kept_taxa/${{gene}}_kept_taxa.txt"
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
        nt="output/after_trimal/outlier_detection/realignment/{gene}.fa",
        aa_msa="output/after_trimal/outlier_detection/realignment/{gene}_aligned.faa",
        pal2nal_path=workflow.basedir+"/utils/extract-buscos/pal2nal.py"
    output:
        nt_aln="output/after_trimal/outlier_detection/realignment/{gene}_aligned.fa"
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
        d=directory("output/genes_to_concat/"),
        okf="output/genes_to_concat/OK"
    shell:
        """
        mkdir -p {output.d}
        cp {input} {output.d}
        touch {output.okf}
        """


rule concatenate:
    input:
        okf="output/genes_to_concat/OK"
    output:
        data="output/supermatrix.phy",
        part="output/supermatrix.nex"
    conda:
        "../envs/phylo_scripts_python.yaml"
    shell:
        """
        rm {input.okf}
        utils/phylo_scripts/concat-aln $(dirname {input.okf}) supermatrix DNA
        mv supermatrix.{{phy,nex}} output/
        """
