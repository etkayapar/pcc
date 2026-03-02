rule init_before_trimal:
    input:
        nt="input/genewise_fastas/{gene}.fa",
        aa="input/genewise_fastas/{gene}.faa"
    output:
        out_nt="output/before_trimal/genewise_fastas/{gene}.fa",
        out_aa="output/before_trimal/genewise_fastas/{gene}.faa",
    group: "first_pass_align"
    shell:
        """
        ln -sr {input.nt} {output.out_nt}
        ln -sr {input.aa} {output.out_aa}
        """

rule align_aa:
    input:
        "output/before_trimal/genewise_fastas/{gene}.faa"
    output:
        "output/before_trimal/mafft_output/{gene}_aligned.faa"
    threads: 4
    group: "first_pass_align"
    conda:
        "../envs/mafft.yaml"
    shell:
        "linsi --thread {threads} --threadit 0 {input} > {output}"

rule backtranslate:
    input:
        nt="output/before_trimal/genewise_fastas/{gene}.fa",
        aa_msa="output/before_trimal/mafft_output/{gene}_aligned.faa",
        pal2nal_path=workflow.basedir+"/utils/extract-buscos/pal2nal.py"
    output:
        "output/before_trimal/mafft_output/{gene}_aligned.fa",
        "output/before_trimal/mafft_output/{gene}.fa",
    log:
        workflow.basedir+"/logs/before_trimal/backtranslate/{gene}_backtranslate.log"
    group: "first_pass_align"
    conda:
        "../envs/phylo_scripts_python.yaml"
    shell:
        """
        ln -r -s {input.nt} output/before_trimal/mafft_output/{wildcards.gene}.fa
        cd output/before_trimal/mafft_output
        (python3 {input.pal2nal_path} {wildcards.gene}) 2> {log}
        """

rule clean_all_gap_seqs:
    input:
        "output/before_trimal/mafft_output/{gene}_aligned.fa"
    output:
        "output/before_trimal/gene_tree_input/{gene}.fa"
    group: "first_pass_align"
    shell:
        "utils/phylo_scripts/cleanAllGaps {input} > {output}"


rule run_treeshrink_before_trimal:
    input:
        treefile="output/before_trimal/outlier_detection/all_genes.treefile",
        gene_names="output/before_trimal/outlier_detection/all_genes_names.txt"
    output:
        treeshrink_output="output/before_trimal/outlier_detection/saved_genes_removed_taxa/output.treefile",
        treeshrink_dir=directory("output/before_trimal/outlier_detection/saved_genes_removed_taxa"),
    log:
        workflow.basedir+"/logs/before_trimal/run_treeshrink.log"
    conda:
        "../envs/treeshrink.yaml"
    params:
        taxa_threshold=config["params"]["detect_outliers"]["taxa_threshold"],
        pipeline_stage="before_trimal",
        treeshrink_mode=config["params"]["detect_outliers"]["treeshrink_mode"],
        long_branch_threshold=config["params"]["detect_outliers"]["long_branch_threshold"]
    script:
        "../utils/detect_outliers_treeshrink.py"

rule detect_outliers_before_trimal:
    input:
        treefile="output/before_trimal/outlier_detection/all_genes.treefile",
        gene_names="output/before_trimal/outlier_detection/all_genes_names.txt",
        treeshrink_output=rules.run_treeshrink_before_trimal.output.treeshrink_output
    output:
        outlier_genes_list="output/before_trimal/outlier_detection/outlier_genes.txt",
    log:
        workflow.basedir+"/logs/before_trimal/detect_outliers.log"
    conda:
        "../envs/detect_outliers.yaml"
    params:
        taxa_threshold=config["params"]["detect_outliers"]["taxa_threshold"],
        pipeline_stage="before_trimal",
        treeshrink_mode=config["params"]["detect_outliers"]["treeshrink_mode"],
        long_branch_threshold=config["params"]["detect_outliers"]["long_branch_threshold"]
    script:
        "../utils/detect_outliers.py"

        
checkpoint process_outliers_before_trimal:
    input:
        fastas=expand("output/before_trimal/gene_tree_input/{gene}.fa", gene=genes),
        treeshrink_dir=rules.run_treeshrink_before_trimal.output.treeshrink_dir,
        remove_taxa_path="utils/phylo_scripts/remove_taxa.awk",
        outlier_genes="output/before_trimal/outlier_detection/outlier_genes.txt"
    output:
        genelist="output/before_trimal/outlier_detection/final_output/genelist.txt",
        d=directory("output/before_trimal/outlier_detection/final_output")
    log:
        workflow.basedir+"/logs/before_trimal/process_outliers.log"
    shell:
        r"""
        set +o pipefail
        (for fasta in {input.fastas}; do
            gene=$(basename $fasta .fa)
            removed_taxa_path="output/before_trimal/outlier_detection/saved_genes_removed_taxa/${{gene}}_removed_taxa.txt"
            if [[ -f $removed_taxa_path ]]
            then
            {input.remove_taxa_path} -v taxafile=${{removed_taxa_path}} $fasta > {output.d}/${{gene}}.fa
            echo ${{gene}} >> {output.genelist}
            elif grep -qw ${{gene}} {input.outlier_genes}
            then
            true
            else
            ln -sr $fasta {output.d}
            echo ${{gene}} >> {output.genelist}
            fi
        done) 2>{log}
        """

rule run_trimal:
    input:
        aa_aln="output/before_trimal/outlier_detection/realignment/{gene}_aligned.faa",
        pal2nal_path=workflow.basedir+"/utils/extract-buscos/pal2nal.py"
    output:
        trimal_cols="output/before_trimal/outlier_detection/realignment/{gene}.trimal",
        nt_aln="output/before_trimal/outlier_detection/realignment/{gene}_aligned.fa"
    conda:
        "../envs/mafft.yaml"
    group: "realign_before_trimal"
    shell:
        """
        trimal -in {input.aa_aln} -out /dev/null -automated1 -colnumbering > {output.trimal_cols}
        cd $(dirname {input.aa_aln})
        {input.pal2nal_path} {wildcards.gene} 
        """
