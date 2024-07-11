include: "rules/00-common.smk"
include: "rules/01-before_trimal.smk"
include: "rules/02-after_trimal.smk"
include: "rules/03-repeated_rules.smk"

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
