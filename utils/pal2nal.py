#! /usr/bin/env python3

import re
import os
import glob
import sys
from functools import reduce
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def parse_trimal(fname):
    '''
    Parse colnumbering notation from TrimAl results
    of individual alignments.
    '''

    import itertools
    def ranges(i):
        for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
            b = list(b)
            yield b[0][1], b[-1][1]

    with open(fname) as trimal:
        trimal_result = trimal.readlines()
        trimal_result = ", ".join([x.strip() for x in trimal_result]).split(", ")
        trimal_result[0] = trimal_result[0].replace("#ColumnsMap\t", "")
        trimal_result = [int(x) + 1 for x in trimal_result if x != '']
        flanks = list(ranges(trimal_result))
    
    try:
        nucl_flanks = [[(x[0]-1)*3, x[1]*3] for x in flanks]
    except TypeError:
        print(f"Error with file: {fname}")

    return nucl_flanks

def parse_gblocks(fname):
    
    with open(fname) as gblocks:
        for l in gblocks:
            if l.startswith("Flanks:"):
                this_flank = l.strip().split(":")[1]
                if this_flank == "":
                    return None
                this_flank = eval(this_flank.replace("  ", ","))
                if any(isinstance(x, list) for x in this_flank):
                    flanks = this_flank
                else:
                    flanks = [this_flank]
    try:
        nucl_flanks = [[(x[0]-1)*3, x[1]*3] for x in flanks]
    except TypeError:
        print(f"Error with file: {fname}")

    return nucl_flanks

def parse_guidance(fname, threshold=0.93):

    '''
    reads the columnwise score file from Guidance2 and discards low confidence
    columns that are below the <threshold> and returns contiguous high
    confidence blocks as list of ranges as ready to use python slicing information
    in nucleotide units.
    '''
    
    flanks = []
    with open(fname) as f:
        tmp_flanks = []
        next(f)
        l = f.readline()
        while not l.startswith("#"):
            score = l.strip().split()
            score[0] = int(score[0])
            ##print(l)
            if float(score[1]) < threshold:
                if len(tmp_flanks) > 0:
                    flanks.append([tmp_flanks[0][0],tmp_flanks[-1][0]])
                    tmp_flanks = []
                l = f.readline()
                continue
            tmp_flanks.append(score)
            l = f.readline()
        if len(tmp_flanks) > 0:
            flanks.append([tmp_flanks[0][0], tmp_flanks[-1][0]])

    
    try:
        nucl_flanks = [[(x[0]-1)*3, x[1]*3] for x in flanks]
    except TypeError as e:
        print(f"Error with file: {fname}")
        print(e)

    return nucl_flanks

def prot_to_codon(prot_alignment, codon_seqs, gap_char="-", codon_table="1"):
    """
    using a protein alignment and unaligned codon sequences,
    return the corresponding nucleotide alignment. This does
    not check for possible intricate problems like frameshifts.
    """


    l = len(prot_alignment)
    if l != len(codon_seqs):
        return

    nucl_records = []

    for sp_idx in range(l):
        this_desc = prot_alignment[sp_idx].description
        this_prot = str(prot_alignment[sp_idx,:].seq)
        this_prot_n = this_prot.replace(gap_char, "")

        this_nucl = codon_seqs[this_desc].seq
        this_nucl_tra = this_nucl.upper().translate(table=codon_table)
        this_nucl = str(this_nucl)

        if this_nucl_tra != this_prot_n:
            ### Print an error/warning
            print(f"[ERROR]: NT and AA seqs are different for the species {this_desc}", file=sys.stderr)
            return

        c = len(this_prot)
        this_nucl_align = ""
        npos = 0

        for aa in this_prot:
            if aa == "-":
                this_nucl_align += "-"*3
                continue
            this_nucl_align += this_nucl[npos:npos+3]
            npos += 3

        this_record = SeqRecord(Seq(this_nucl_align), id=this_desc.split("|")[0],
                                description="")

        nucl_records.append(this_record)

    return MultipleSeqAlignment(nucl_records)

def parse_alignment(gene, codon_table="1"):
    """
    Parse the prot alignment, parse the unaligned codon fasta
    rebuild the codon alignment from them. Then, parse the
    Gblocks output (coordinates of conserved regions) and
    subset the codon alignment accordingly.
    """
    prot_alignment_fname = gene+"_aligned.faa"
    nucl_unaligned_fname = gene+".fa"
    gblocks_output_fname = gene+"_aligned.faa-gb.htm"
    guidance_output_fname= gene+"_guidance_col_col.scr"
    trimal_output_fname   = gene+".trimal"

    if os.path.exists(gblocks_output_fname):
        nucl_flanks = parse_gblocks(gblocks_output_fname)
        if not nucl_flanks: ## NO remaining columns after filtering
            return None

    if os.path.exists(guidance_output_fname):
        nucl_flanks = parse_guidance(guidance_output_fname)
        if not nucl_flanks: ## NO remaining columns after filtering
            return None

    if os.path.exists(trimal_output_fname):
        nucl_flanks = parse_trimal(trimal_output_fname)
        if not nucl_flanks: ## NO remaining columns after filtering
            return None

    prot_alignment = AlignIO.read(prot_alignment_fname, "fasta")
    nucl_seqs = SeqIO.to_dict(SeqIO.parse(nucl_unaligned_fname, "fasta"), key_function = lambda rec: rec.description)
    codon_alignment = prot_to_codon(prot_alignment, nucl_seqs, codon_table=codon_table)

    if os.path.exists(gblocks_output_fname) or os.path.exists(guidance_output_fname) or os.path.exists(trimal_output_fname):
        try:
            clean_alignment = reduce(lambda a,b: a + b, [codon_alignment[:,x[0]:x[1]] for x in nucl_flanks])
        except TypeError:
            exit(1)
    else:
        clean_alignment = codon_alignment
    
    clean_alignment.sort()

    return clean_alignment

def usage():
    msg = """Bactranslate amino-acid alignments to nucleotide alignments with the optional incorporation of alignment
trimming information
---------------------------------------------------
Usage:
pal2nal.py <GENE_BASENAME> [<NCBI_CODON_TABLE>]

ARGUMENTS:
    GENE_BASENAME             Filename prefix for the gene name without the .fa extension. Hypothetically if the
                            value for this argument is "gene", then the following files should be present in the
                            directory where this script is called from: 1) gene.fa (unaligned NT seqs),
                            2) gene_aligned.faa (AA MSA), and 3) gene.trimal (comma separated columns to retain,
                            trimal -colnumbering output) . The last file is optional in case you want to backtranslate
                            while incorporating alignment trimming information. After a successful run, the program
                            outputs a new file to the same directory called "gene_aligned.fa".

    NCBI_CODON_TABLE <INT>    Integer code for the desired NCBI genetic code table [default: 1 (The standard code)]
                            use 5 for invertebrate mitochondrial code.
          """

    print(msg, file=sys.stderr)
    exit(1)

def main():
    try:
        gene=sys.argv[1]
    except IndexError:
        usage()
    if gene=="-h":
        usage()
    try:
        codon_table=sys.argv[2]
    except IndexError:
        codon_table=1

    nucl_aln = parse_alignment(gene, codon_table)

    AlignIO.write(nucl_aln, gene+"_aligned.fa", "fasta-2line")

if __name__ =="__main__":
    main()
