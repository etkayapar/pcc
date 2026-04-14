containerized: "docker://ghcr.io/etkayapar/pcc:0.1.0"

include: "rules/00-common.smk"
include: "rules/01-before_trimal.smk"
include: "rules/02-after_trimal.smk"
include: "rules/03-repeated_rules.smk"

localrules: init_before_trimal,init_after_trimal,run_trimal,backtranslate,backtranslate_final,clean_all_gap_seqs,unalign_outliers_before_trimal,unalign_outliers_after_trimal

## Target rules for a disjunct, three-step workflow
rule first_pass:
    input:
        "output/before_trimal/outlier_detection/outlier_genes.txt"

rule second_pass:
    input:
        "output/after_trimal/outlier_detection/outlier_genes.txt"

rule conclude:
    input:
        "output/supermatrix.phy"

rule conclude_with_gene_trees:
    input:
        "output/supermatrix.phy",
        expand(
            "output/final_gene_trees/{gene}/{gene}.treefile",
            gene=get_filtered_genes_final
        ),
        
