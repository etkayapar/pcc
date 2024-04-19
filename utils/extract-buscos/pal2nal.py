#! /usr/bin/env python3

import find_buscos
import sys

gene=sys.argv[1]
try:
    codon_table=sys.argv[2]
except IndexError:
    codon_table=1

nucl_aln = find_buscos.parse_alignment(gene, codon_table)

find_buscos.AlignIO.write(nucl_aln, gene+"_aligned.fa", "fasta-2line")
