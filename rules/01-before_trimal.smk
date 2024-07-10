rule init_before_trimal:
    input:
        nt="genewise_fastas/{gene}.fa",
        aa="genewise_fastas/{gene}.faa"
    output:
        out_nt="before_trimal/genewise_fastas/{gene}.fa",
        out_aa="before_trimal/genewise_fastas/{gene}.faa"
    shell:
        """
        ln -sr {input.nt} {output.out_nt}
        ln -sr {input.aa} {output.out_aa}
        """

rule align_aa:
    input:
        "before_trimal/genewise_fastas/{gene}.faa"
    output:
        "before_trimal/mafft_output/{gene}_aligned.faa"
    threads: 4
    conda:
        "../envs/mafft.yaml"
    shell:
        "linsi --thread {threads} {input} > {output}"

rule backtranslate:
    input:
        nt="before_trimal/genewise_fastas/{gene}.fa",
        aa_msa="before_trimal/mafft_output/{gene}_aligned.faa",
        pal2nal_path=workflow.basedir+"/utils/extract-buscos/pal2nal.py"
    output:
        "before_trimal/mafft_output/{gene}_aligned.fa"
    log:
        workflow.basedir+"/logs/before_trimal/backtranslate/{gene}_backtranslate.log"
    conda:
        "../envs/phylo_scripts_python.yaml"
    shell:
        """
        ln -r -s {input.nt} before_trimal/mafft_output/{wildcards.gene}.fa
        cd before_trimal/mafft_output
        (python3 {input.pal2nal_path} {wildcards.gene}) 2> {log}
        """

rule clean_all_gap_seqs:
    input:
        "before_trimal/mafft_output/{gene}_aligned.fa"
    output:
        "before_trimal/gene_tree_input/{gene}.fa"
    shell:
        "utils/phylo_scripts/cleanAllGaps {input} > {output}"


rule detect_outliers_before_trimal:
    input:
        treefile="before_trimal/outlier_detection/all_genes.treefile",
        gene_names="before_trimal/outlier_detection/all_genes_names.txt"
    output:
        saved_genes_plot  ="before_trimal/outlier_detection/saved_genes.pdf",
        outlier_genes_plot="before_trimal/outlier_detection/outlier_genes.pdf",
        outlier_genes_list="before_trimal/outlier_detection/outlier_genes.txt",
        kept_taxa_dir=directory("before_trimal/outlier_detection/saved_genes_kept_taxa")
    conda:
        "../envs/detect_outliers.yaml"
    params:
        long_branch_threshold=0.1,
        taxa_threshold=0.7
    shell:
        f"""
        Rscript utils/gene_trees.R {{params.long_branch_threshold}} {{params.taxa_threshold}} {{input.treefile}} {{input.gene_names}} before_trimal/outlier_detection > gt.log 2>gt.err
        """

checkpoint process_outliers_before_trimal:
    input:
        aln_dir="before_trimal/gene_tree_input",
        keep_taxa_path="utils/phylo_scripts/keep_taxa.awk"
    output:
        genelist="before_trimal/outlier_detection/final_output/genelist.txt",
        d=directory("before_trimal/outlier_detection/final_output")
        #"{stage}/outlier_detection/final_output/{gene}.fa"
    params:
        #kept_taxa_path="before_trimal/outlier_detection/saved_genes_kept_taxa/{{gene}}_kept_taxa.txt",
        outlier_genes_path="before_trimal/outlier_detection/outlier_genes.txt"
    log:
        workflow.basedir+"/logs/before_trimal/process_outliers.log"
    shell:
        """
        set +o pipefail
        (for gene in `ls {input.aln_dir}/ | grep -E "*.fa$" | cut -d. -f 1`
        do
            kept_taxa_path="before_trimal/outlier_detection/saved_genes_kept_taxa/${{gene}}_kept.taxa.txt"
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

rule run_trimal:
    input:
        aa_aln="before_trimal/outlier_detection/realignment/{gene}_aligned.faa",
        pal2nal_path=workflow.basedir+"/utils/extract-buscos/pal2nal.py"
    output:
        trimal_cols="before_trimal/outlier_detection/realignment/{gene}.trimal",
        nt_aln="before_trimal/outlier_detection/realignment/{gene}_aligned.fa"
    conda:
        "../envs/mafft.yaml"
    shell:
        """
        trimal -in {input.aa_aln} -out /dev/null -automated1 -colnumbering > {output.trimal_cols}
        cd $(dirname {input.aa_aln})
        {input.pal2nal_path} {wildcards.gene} 
        """
