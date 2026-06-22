STAGES=["before_trimal", "after_trimal"]
for stage in STAGES:
    if config["params"]["intermediate_tree_method"] == "iqtree":
        rule:
            name: f"infer_gene_trees_{stage}_iqtree"
            input:
                f"output/{stage}/gene_tree_input/{{gene}}.fa"
            output:
                treefile=f"output/{stage}/gene_trees/{{gene}}/{{gene}}.treefile",
            threads: 4
            conda:
                "../envs/iqtree.yaml"
            params:
                prefix=f"output/{stage}/gene_trees/{{gene}}/{{gene}}"
            resources:
                runtime="2d"
            shell:
                """
                iqtree2 -s {input} \
                    -m MFP -mset GTR -mrate I+R \
                    -T {threads} --prefix {params.prefix} -st DNA --keep-ident
                """

    elif config["params"]["intermediate_tree_method"] == "fasttree":
        rule:
            name: f"infer_gene_trees_{stage}_fasttree"
            input:
                f"output/{stage}/gene_tree_input/{{gene}}.fa"
            output:
                treefile=f"output/{stage}/gene_trees/{{gene}}/{{gene}}.treefile",
            threads: 1
            conda:
                "../envs/fasttree.yaml"
            resources:
                runtime="2h"
            shell:
                """
                fasttree -gtr -gamma -nt < {input} > {output.treefile}
                """

    else:
        raise ValueError(
                "Set params > intermediate_tree_method in the config.yaml to "
                "either 'fasttree' or 'iqtree'"
            )

    rule:
        name: f"collect_gene_trees_{stage}"
        input:
            branch(stage=="before_trimal",
                   then=expand(f"output/{stage}/gene_trees/{{gene}}/{{gene}}.treefile", gene=genes),
                   otherwise=get_gene_list_to_infer_tree_after)
        output:
            trees=f"output/{stage}/outlier_detection/all_genes.treefile",
            gene_names=f"output/{stage}/outlier_detection/all_genes_names.txt"
        shell:
            """
            find {input} | cut -d '/' -f4 > {output.gene_names}
            cat {input} > {output.trees}
            """
    rule:
        name: f"unalign_outliers_{stage}"
        input:
            nt=f"output/{stage}/outlier_detection/final_output/{{gene}}.fa",
            unaln_scr_path=workflow.source_path("../utils/phylo_scripts/unalignFasta.awk"),
            tra_path=workflow.source_path("../utils/phylo_scripts/translate_stdin.py"),
            gaps_scr_path=workflow.source_path("../utils/phylo_scripts/cleanAllGaps")
        output:
            nt_unaligned=f"output/{stage}/outlier_detection/realignment/{{gene}}.fa",
            aa_unaligned=f"output/{stage}/outlier_detection/realignment/{{gene}}.faa"
        wildcard_constraints:
            gene="[A-Za-z0-9]+"
        conda:
            "../envs/phylo_scripts_python.yaml"
        params:
            maxgap_pct=config["params"]["general"]["maxgap_pct"]
        shell:
            """
            awk -f {input.gaps_scr_path} -v p={params.maxgap_pct} {input.nt} | awk -f {input.unaln_scr_path}  > {output.nt_unaligned}
            cat {output.nt_unaligned} | python3 {input.tra_path}  > {output.aa_unaligned}
            """
    rule:
        name: f"realign_outliers_{stage}"
        input:
            aa=f"output/{stage}/outlier_detection/realignment/{{gene}}.faa"
        output:
            aa_aln=ensure(f"output/{stage}/outlier_detection/realignment/{{gene}}_aligned.faa", non_empty=True)
        conda:
            "../envs/mafft.yaml"
        threads: 4
        params:
            aligner=get_aln_params,
            mafft_tmpdir=config["params"]["align_aa"]["mafft_tmpdir"],
            stage=stage
        shell:
            """
            mafft_tmpdir={params.mafft_tmpdir}
            if [[ -n ${{mafft_tmpdir}} ]]
            then
                mafft_tmpdir=${{mafft_tmpdir}}/realign_outliers_{params.stage}/{wildcards.gene}
                mkdir -p $mafft_tmpdir
            fi

            MAFFT_TMPDIR=$mafft_tmpdir {params.aligner} --thread {threads} {input.aa} > {output.aa_aln}

            if [[ -n ${{mafft_tmpdir}} ]]
            then
                rm -rf ${{mafft_tmpdir}}
            fi
            """
