include: "common_utils.smk"
include: "rules/00-common.smk"
include: "rules/01-before_trimal.smk"

## Target rules for a disjunct, three-step workflow
rule first_pass:
    input:
        "before_trimal/outlier_detection/outlier_genes.txt"

rule second_pass:
    input:
        "after_trimal/outlier_detection/outlier_genes.txt"

rule conclude:
    input:
        "supermatrix.phy"
