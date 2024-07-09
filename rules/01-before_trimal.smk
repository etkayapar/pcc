stage = "before_trimal"

rule first_pass:
    input:
        "before_trimal/outlier_detection/outlier_genes.txt"

rule second_pass:
    input:
        "after_trimal/outlier_detection/outlier_genes.txt"

rule conclude:
    input:
        "supermatrix.phy"

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

rule:
    name: f"infer_gene_trees_{stage}"
    input:
        f"{stage}/gene_tree_input/{{gene}}.fa"
    output:
        treefile=f"{stage}/gene_trees/{{gene}}/{{gene}}.treefile",
        treedir=directory(f"{stage}/gene_trees/{{gene}}")
    threads: 4
    conda:
        "../envs/iqtree.yaml"
    shell:
        """
        iqtree2 -s {input} \
                -m MFP -mset GTR -mrate I+R \
                -T {threads} --prefix {wildcards.gene} -st DNA --keep-ident
        mv {wildcards.gene}.* {output.treedir}/
        """

rule:
    name: f"collect_gene_trees_{stage}"
    input:
        expand(f"{stage}/gene_trees/{{gene}}/{{gene}}.treefile", gene=genes)
        #get_gene_list
    output:
        trees=f"{stage}/outlier_detection/all_genes.treefile",
        gene_names=f"{stage}/outlier_detection/all_genes_names.txt"
    shell:
        """
        find {input} | cut -d '/' -f3 > {output.gene_names}
        cat {input} > {output.trees}
        """

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
        aln_dir=f"{stage}/gene_tree_input",
        keep_taxa_path="utils/phylo_scripts/keep_taxa.awk"
    output:
        genelist=f"{stage}/outlier_detection/final_output/genelist.txt",
        d=directory(f"{stage}/outlier_detection/final_output")
        #"{stage}/outlier_detection/final_output/{gene}.fa"
    params:
        #kept_taxa_path=f"{stage}/outlier_detection/saved_genes_kept_taxa/{{gene}}_kept_taxa.txt",
        outlier_genes_path=f"{stage}/outlier_detection/outlier_genes.txt"
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
# checkpoint process_outliers_before_trimal:
#     input:
#         aln_path=f"{stage}/gene_tree_input/{{gene}}.fa",
#         keep_taxa_path="utils/phylo_scripts/keep_taxa.awk"
#     output:
#         genelist=f"{stage}/outlier_detection/final_output/genelist.txt",
#         d=directory(f"{stage}/outlier_detection/final_output")
#         #"{stage}/outlier_detection/final_output/{gene}.fa"
#     params:
#         kept_taxa_path=f"{stage}/outlier_detection/saved_genes_kept_taxa/{{gene}}_kept_taxa.txt",
#         outlier_genes_path=f"{stage}/outlier_detection/outlier_genes.txt"
#     shell:
#         """
#         set +o pipefail
#         if [[ -f {params.kept_taxa_path} ]]
#         then
#         {input.keep_taxa_path} -v taxafile={params.kept_taxa_path} {input.aln_path} > {output.d}/{wildcards.gene}.fa
#         echo {wildcards.gene} >> {output.genelist}
#         elif grep -qw {wildcards.gene} {params.outlier_genes_path}
#         then
#         true
#         else
#         ln -sr {input.aln_path} {output.d}
#         echo {wildcards.gene} >> {output.genelist}
#         fi
#         """

rule:
    name: f"unalign_outliers_{stage}"
    input:
        nt=f"{stage}/outlier_detection/final_output/{{gene}}.fa",
        scr_path="utils/phylo_scripts/unalignFasta.awk",
        tra_path="utils/phylo_scripts/translate_stdin.py"
    output:
        nt_unaligned=f"{stage}/outlier_detection/realignment/{{gene}}.fa",
        aa_unaligned=f"{stage}/outlier_detection/realignment/{{gene}}.faa"
    wildcard_constraints:
        gene="[A-Za-z0-9]+"
    conda:
        "../envs/phylo_scripts_python.yaml"
    shell:
        """
        {input.scr_path} {input.nt} > {output.nt_unaligned}
        cat {output.nt_unaligned} | {input.tra_path}  > {output.aa_unaligned}
        """

rule:
    name: f"realign_outliers_{stage}"
    input:
        aa=f"{stage}/outlier_detection/realignment/{{gene}}.faa"
    output:
        aa_aln=f"{stage}/outlier_detection/realignment/{{gene}}_aligned.faa"
    conda:
        "../envs/mafft.yaml"
    threads: 4
    shell:
        """
        linsi --thread {threads} {input.aa} > {output.aa_aln}
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

stage="after_trimal"
rule init_after_trimal:
    input:
        "before_trimal/outlier_detection/realignment/{gene}_aligned.fa"
    output:
        "after_trimal/gene_tree_input/{gene}.fa"
    shell:
        """
        utils/phylo_scripts/cleanAllGaps {input} > {output}
        """

rule:
    name: f"infer_gene_trees_{stage}"
    input:
        f"{stage}/gene_tree_input/{{gene}}.fa"
    output:
        treefile=f"{stage}/gene_trees/{{gene}}/{{gene}}.treefile",
        treedir=directory(f"{stage}/gene_trees/{{gene}}")
    threads: 4
    conda:
        "../envs/iqtree.yaml"
    shell:
        """
        iqtree2 -s {input} \
                -m MFP -mset GTR -mrate I+R \
                -T {threads} --prefix {wildcards.gene} -st DNA --keep-ident
        mv {wildcards.gene}.* {output.treedir}/
        """

rule:
    name: f"collect_gene_trees_{stage}"
    input:
        get_gene_list_to_infer_tree
        ##input func on checkpoint
        ##expand(f"{stage}/gene_trees/{{gene}}/{{gene}}.treefile", gene=genes)
    output:
        trees=f"{stage}/outlier_detection/all_genes.treefile",
        gene_names=f"{stage}/outlier_detection/all_genes_names.txt"
    shell:
        """
        find {input} | cut -d '/' -f3 > {output.gene_names}
        cat {input} > {output.trees}
        """

rule:
    name: f"detect_outliers_{stage}"
    input:
        treefile=f"{stage}/outlier_detection/all_genes.treefile",
        gene_names=f"{stage}/outlier_detection/all_genes_names.txt"
    output:
        saved_genes_plot  =f"{stage}/outlier_detection/saved_genes.pdf",
        outlier_genes_plot=f"{stage}/outlier_detection/outlier_genes.pdf",
        outlier_genes_list=f"{stage}/outlier_detection/outlier_genes.txt",
        kept_taxa_dir=directory(f"{stage}/outlier_detection/saved_genes_kept_taxa")
    conda:
        "../envs/detect_outliers.yaml"
    params:
        long_branch_threshold=0.1,
        taxa_threshold=0.7
    shell:
        f"""
        Rscript utils/gene_trees.R {{params.long_branch_threshold}} {{params.taxa_threshold}} {{input.treefile}} {{input.gene_names}} {stage}/outlier_detection > gt.log 2>gt.err
        """

checkpoint process_outliers_after_trimal:
    input:
        aln_dir=f"{stage}/gene_tree_input",
        keep_taxa_path="utils/phylo_scripts/keep_taxa.awk"
    output:
        genelist=f"{stage}/outlier_detection/final_output/genelist.txt",
        d=directory(f"{stage}/outlier_detection/final_output")
        #"{stage}/outlier_detection/final_output/{gene}.fa"
    params:
        #kept_taxa_path=f"{stage}/outlier_detection/saved_genes_kept_taxa/{{gene}}_kept_taxa.txt",
        outlier_genes_path=f"{stage}/outlier_detection/outlier_genes.txt"
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

rule:
    name: f"unalign_outliers_{stage}"
    input:
        nt=f"{stage}/outlier_detection/final_output/{{gene}}.fa",
        scr_path="utils/phylo_scripts/unalignFasta.awk",
        tra_path="utils/phylo_scripts/translate_stdin.py"
    output:
        nt_unaligned=f"{stage}/outlier_detection/realignment/{{gene}}.fa",
        aa_unaligned=f"{stage}/outlier_detection/realignment/{{gene}}.faa"
    conda:
        "../envs/phylo_scripts_python.yaml"
    shell:
        """
        {input.scr_path} {input.nt} > {output.nt_unaligned}
        cat {output.nt_unaligned} | {input.tra_path}  > {output.aa_unaligned}
        """

rule:
    name: f"realign_outliers_{stage}"
    input:
        aa=f"{stage}/outlier_detection/realignment/{{gene}}.faa"
    output:
        aa_aln=f"{stage}/outlier_detection/realignment/{{gene}}_aligned.faa"
    conda:
        "../envs/mafft.yaml"
    threads: 4
    shell:
        """
        linsi --thread {threads} {input.aa} > {output.aa_aln}
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
             
