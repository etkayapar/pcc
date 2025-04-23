rule init_after_trimal:
    input:
        "output/before_trimal/outlier_detection/realignment/{gene}_aligned.fa"
    output:
        "output/after_trimal/gene_tree_input/{gene}.fa"
    shell:
        """
        utils/phylo_scripts/cleanAllGaps {input} > {output}
        """

rule run_treeshrink_after_trimal:
    input:
        treefile="output/after_trimal/outlier_detection/all_genes.treefile",
        gene_names="output/after_trimal/outlier_detection/all_genes_names.txt"
    output:
        treeshrink_output="output/after_trimal/outlier_detection/saved_genes_removed_taxa/output.treefile"
    conda:
        "../envs/treeshrink.yaml"
    params:
        taxa_threshold=config["params"]["detect_outliers"]["taxa_threshold"],
        pipeline_stage="after_trimal",
        treeshrink_mode=config["params"]["detect_outliers"]["treeshrink_mode"],
        long_branch_threshold=config["params"]["detect_outliers"]["long_branch_threshold"]
    script:
        "../utils/detect_outliers_treeshrink.py"

rule detect_outliers_after_trimal:
    input:
        treefile="output/after_trimal/outlier_detection/all_genes.treefile",
        gene_names="output/after_trimal/outlier_detection/all_genes_names.txt",
        treeshrink_output=rules.run_treeshrink_after_trimal.output.treeshrink_output
    output:
        outlier_genes_list="output/after_trimal/outlier_detection/outlier_genes.txt",
    conda:
        "../envs/detect_outliers.yaml"
    params:
        taxa_threshold=config["params"]["detect_outliers"]["taxa_threshold"],
        pipeline_stage="after_trimal",
        treeshrink_mode=config["params"]["detect_outliers"]["treeshrink_mode"],
        long_branch_threshold=config["params"]["detect_outliers"]["long_branch_threshold"]
    script:
        "../utils/detect_outliers.py"

checkpoint process_outliers_after_trimal:
    input:
        aln_dir="output/after_trimal/gene_tree_input",
        remove_taxa_path="utils/phylo_scripts/remove_taxa.awk"
    output:
        genelist="output/after_trimal/outlier_detection/final_output/genelist.txt",
        d=directory("output/after_trimal/outlier_detection/final_output")
    params:
        outlier_genes_path="output/after_trimal/outlier_detection/outlier_genes.txt"
    log:
        workflow.basedir+"/logs/after_trimal/process_outliers.log"
    shell:
        """
        set +o pipefail
        (for gene in `ls {input.aln_dir}/ | grep -E ".*\.fa$" | cut -d. -f 1`
        do
            removed_taxa_path="output/after_trimal/outlier_detection/saved_genes_removed_taxa/${{gene}}_removed_taxa.txt"
            if [[ -f $removed_taxa_path ]]
            then
            {input.remove_taxa_path} -v taxafile=${{removed_taxa_path}} {input.aln_dir}/${{gene}}.fa > {output.d}/${{gene}}.fa
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
        genes=get_gene_list_to_concatenate,
        get_aln_len="utils/phylo_scripts/get_aln_len.awk"
    output:
        d=directory("output/genes_to_concat/"),
        okf="output/genes_to_concat/OK",
        shortgenes="output/too_short_genes.txt"
    params:
        min_aln_len=config["params"]["concatenate"]["min_aln_len"]
    shell:
        """
        mkdir -p {output.d}
        touch {output.shortgenes}
        for gene in {input.genes}
        do
            aln_len=$({input.get_aln_len} ${{gene}})
            if [ $aln_len -ge {params.min_aln_len} ] 
            then
                cp ${{gene}} {output.d}
            else
                echo ${{gene}} >> {output.shortgenes}
            fi
        done
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
        if [ -z "$( ls $(dirname {input.okf}))" ]
        then
        echo "ERROR: None of the genes passed the minimum alignment length threshold"
        exit 1
        else
        utils/phylo_scripts/concat-aln $(dirname {input.okf}) supermatrix DNA
        mv supermatrix.{{phy,nex}} output/
        fi
        """
