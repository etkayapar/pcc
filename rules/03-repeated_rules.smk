STAGES=["before_trimal", "after_trimal"]
for stage in STAGES:
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
            branch(stage=="before_trimal",
                   then=expand(f"{stage}/gene_trees/{{gene}}/{{gene}}.treefile", gene=genes),
                   otherwise=get_gene_list_to_infer_tree_after)
        output:
            trees=f"{stage}/outlier_detection/all_genes.treefile",
            gene_names=f"{stage}/outlier_detection/all_genes_names.txt"
        shell:
            """
            find {input} | cut -d '/' -f3 > {output.gene_names}
            cat {input} > {output.trees}
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
