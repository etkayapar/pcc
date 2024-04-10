rule align_aa:
    input:
        "genewise_fastas/{gene}.faa"
    output:
        "mafft_output/{gene}_aligned.faa"
    threads: 8
    conda:
        "envs/mafft.yaml"
    shell:
        "linsi --thread {threads} {input} > {output}"

rule backtranslate:
    input:
        nt="genewise_fastas/{gene}.fa",
        aa_msa="mafft_output/{gene}_aligned.faa"
    output:
        "mafft_output/{gene}_aligned.fa"
    shell:
        """
        ln -s {input.nt} mafft_output/{gene}.fa
        cd mafft_output
        python3 utils/extract-buscos/pal2nal.py {gene}
        """
rule clean_all_gap_seqs:
    input:
        "mafft_output/{gene}_aligned.fa"
    output:
        "before_trimal/{gene}.fa"
    shell:
        "utils/phylo_scripts/cleanAllGaps {input} > {output}"

rule infer_gene_trees:
    input:
        "before_trimal/{gene}.fa"
    output:
        "gene_trees/before_trimal/{gene}.treefile"
    threads: 4
    conda:
        "envs/iqtree.yaml"
    shell:
        """
        iqtree2 -s {input} \
                -m MFP -mset GTR -mrate I+R \
                -T {threads} --prefix {gene} -st DNA --keep-ident
        mkdir -p {gene}
        mv {gene}.* {gene}/
        """
        
        
    
