import re
import os
import glob
from functools import reduce
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from collections import Counter
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Nexus import Nexus

deneme_list = ["TRINITY_DN161_c0_g2_i13:2390-28141",
               "TRINITY_DN161_c0_g2_i13:2390-28141"]
BUSCO_DATASET = "run_lepidoptera_odb10"

def detect_busco_mode(startdir):
    species_list = [ os.path.basename(f.path) for f in os.scandir(startdir) if (f.is_dir() and "busco_downloads" not in f.path and os.path.basename(f.path)[0] != ".")] 
    retval = {}
    for species in species_list:
        try:
            log_fname = [x for x in glob.glob(startdir+"/"+species+"/short_summary*.txt")][0]
        except IndexError:
            ## The 'species' folder is not a busco_output.
            continue
        with open(log_fname) as logF:
            for l in logF:
                if l.startswith("# BUSCO was run in mode:"):
                    tmpval = l.strip().split(":")[1].lstrip()
                    if "_" in tmpval:
                        tmpval = tmpval.split("_")[1]
                    retval[species] = tmpval
    return retval


def parse_table_transcriptome(fname):
    busco_dict = {}
    with open(fname) as table:
        for l in table:
            if l.startswith("#"):
                continue
            l = l.strip().split()
            if l[1] in ("Missing"):
                continue
            busco = l[0]
            # try:
            pos = l[2].split(":")[1].split("-")
            # except IndexError:
            #     print(fname)
            #     print(l)
            #     raise IndexError
            l[2] = l[2].split(":")[0]
            rec_id = "|".join([l[0],l[2]]+pos)
            thisEntry = dict(zip(["status", "sequence", "score", "length", "rec_id"], l[1:5]+[rec_id]))
            try:
                busco_dict[busco].append(thisEntry)
            except KeyError:
                busco_dict[busco] = []
                busco_dict[busco].append(thisEntry)
    return busco_dict
                
def parse_table_genome(fname):
    busco_dict = {}
    with open(fname) as table:
        for l in table:
            if l.startswith("#"):
                continue
            l = l.strip().split()
            if l[1] in ("Missing", "Duplicated"):
                continue
            busco = l[0]
            l[2] = l[2].split(":")[0]
            rec_id = "|".join([l[0]] + l[2:5])
            thisEntry = dict(zip(["status", "sequence", "score", "length", "rec_id"], l[1:3]+l[6:8]+[rec_id]))
            busco_dict[busco] = thisEntry
    return busco_dict

def parse_table_ids(fname):
    '''
    Parse one of the tables at random so that there will be full list
    BUSCOS analysed in the memory for later use
    '''
    ids = set()
    with open(fname) as table:
        for l in table:
            if l.startswith("#"):
                continue
            l = l.strip().split()
            busco = l[0]
            ids.add(busco)

    ids = list(ids)
    ids.sort()
    return ids

def parse_table(fname, mode):
    cmdtable = {"genome": parse_table_genome,
                "transcriptome": parse_table_transcriptome,
                "ids": parse_table_ids}

    return cmdtable[mode](fname)


def parse_fasta(species, busco_dict, basedir):

    def get_overlap_len(coordinates):

        '''
        return the overlap range btw two features

        <coordinates> should be a tuple of tuples with a structure like below.
        (block_min, block_max), (match['sstart'], match['send'])
        '''
        overlap = (max(coordinates[0][0], coordinates[1][0]), min(coordinates[0][1],coordinates[1][1]))

        return overlap[1] - overlap[0] + 1

    fasta_dict = {}
    dirs = ["initial_results", "rerun_results"]
    suffs = ["*.fna", "*.fa", "*.fasta"]

    suffixes = [ os.path.join(d,suf+".codon.fas") for d in dirs for suf in suffs]

    for suffix in suffixes:
        #fname_pat = "./" + species + suffix
        fname_pat = os.path.join(basedir, species, BUSCO_DATASET, "metaeuk_output", suffix)
        #print(fname_pat)
        try:
            fname = [ x for x in glob.glob(fname_pat)][0]
        except IndexError:
            continue
            #raise FileNotFoundError("There is no fasta file at the path: " + fname_pat)

        for record in SeqIO.parse(fname, "fasta"):
            rec_id = record.id.split("|")
            rec_id[0] = rec_id[0].split("_")[0]
            busco = rec_id[0]
            rec_id = "|".join(rec_id[0:2] + rec_id[6:8])
            try:
                if busco_dict[busco][species] == rec_id:
                    seq = str(record.seq)
                    prot= str(record.seq.translate())
                    thisEntry = dict(zip(["busco", "gene", "seq", "prot"], [busco, rec_id, seq, prot]))
                    fasta_dict[busco] = thisEntry
                elif rec_id.split("|")[1] == busco_dict[busco][species].split("|")[1] and get_overlap_len([[int(x) for x in busco_dict[busco][species].split("|")[2:4]], [int(x) for x in rec_id.split("|")[2:4]]]) > 0:
                    seq = str(record.seq)
                    prot= str(record.seq.translate())
                    thisEntry = dict(zip(["busco", "gene", "seq", "prot"], [busco, rec_id, seq, prot]))
                    fasta_dict[busco] = thisEntry
                    
                elif busco_dict[busco][species] == '':
                    gene = ''
                    seq = ''
                    prot = ''
                    thisEntry = dict(zip(["busco", "gene", "seq", "prot"], [busco, gene, seq, prot]))
                    fasta_dict[busco] = thisEntry
                #else:
                    # print("---------")
                    # print(busco_dict[busco][species])
                    # print(record.id)
                    # print("---------")

            except KeyError:
                continue ## That busco is not present anyways
            except IndexError as e:
                continue
                # print(species)
                # print(fname)
                # raise IndexError(e)
    
    return fasta_dict

def parse_fasta_old(species, busco_dict):
    fasta_dict = {}
    suffixes = ["initial_results/*.fna.codon.fas", "initial_results/*.fa.codon.fas", "rerun_results/*.fna.codon.fas", "rerun_results/*.fa.codon.fas"]
    for suffix in suffixes:
        fname_pat = species + "/" + BUSCO_DATASET + "/metaeuk_output/" + suffix
        try:
            fname = [ x for x in glob.glob(fname_pat)][0]
        except IndexError:
            continue
            #raise FileNotFoundError("There is no fasta file at the path: " + fname_pat)

        for record in SeqIO.parse(fname, "fasta"):
            rec_id = record.id.split("|")
            rec_id[0] = rec_id[0].split("_")[0]
            busco = rec_id[0]
            try:
                if not rec_id[1] in busco_dict[busco][species]:
                    if busco_dict[busco][species] == '':
                        gene = ''
                        seq = ''
                        prot = ''
                        thisEntry = dict(zip(["busco", "gene", "seq", "prot"], [busco, gene, seq, prot]))
                        fasta_dict[busco] = thisEntry
                    continue
            except KeyError:
                continue ## That busco is not present anyways
            except IndexError:
                print(species)
                print(fname)
                exit()
            seq = str(record.seq)
            prot= str(record.seq.translate())
            thisEntry = dict(zip(["busco", "gene", "seq", "prot"], [busco, rec_id[1], seq, prot]))
            ##try:
            fasta_dict[busco] = thisEntry
                ##fasta_dict[busco].append(thisEntry)
            ##except KeyError:
                ##pass ### that busco is not present anwyays
                # fasta_dict[busco] = []
                # fasta_dict[busco].append(thisEntry)

    
    return fasta_dict

def is_same_gene(seqs, delimiter):
    '''
   Checks a single busco_id that is reported as duplicated in the
   full_table.csv to see if multiple entries are actually multi-hits
   or isoforms of the same gene (Happens in case of transcriptome runs)

   busco_id: the id for the busco gene (looks like 1234at213), or the value
             in the first column of full_table.csv
   
   data: the data that stores the information in the full_table.csv

   delimiter: How to parse the sequence id (3rd column) to determine if the 
              duplicates belong to the same gene. (Maybe an actual regex 
              pattern?)
   
   returns: True if duplicates are the same gene
    '''

    if delimiter == "trinity":
        delimiter = "(TRINITY_DN[0-9]+_c[0-9]+)_(g[0-9]+)_(i[0-9]+)"
    gene_pat = re.compile(delimiter)
    try:
        gene_ids = { re.findall(gene_pat, seq_id)[0][1] for seq_id in seqs }
    except IndexError:
        delimiter = "(comp[0-9]+)_(c[0-9]+)_(seq[0-9]+)"
        gene_pat = re.compile(delimiter)
        gene_ids = { re.findall(gene_pat, seq_id)[0][1] for seq_id in seqs }
        print(seqs)

    return len(gene_ids) == 1

def find_best_duplicate(busco_id, full_table):
    '''
    Find and retrieve the best hit among the duplicates.
    Running this function all busco_ids in effect reduces
    each element to a dict (rather than a list of dicts)
    '''
    this_busco = full_table[busco_id]
    this_seqs = [x["sequence"] for x in this_busco]
    if not is_same_gene(this_seqs, "trinity"):
        return None
    scores = [x["score"] for x in this_busco]
    max_idx = scores.index(max(scores))
    return this_busco[max_idx]

def filter_full_table(full_table):
    this_table = full_table.copy()
    for busco_id in list(this_table.keys()):
        filtered = find_best_duplicate(busco_id, this_table)
        if not filtered:
            del(this_table[busco_id])
            continue
        filtered["fastaheader"] = ">" + busco_id + "|" + "_".join(filtered["sequence"].split("_")[0:4])
        this_table[busco_id] = filtered

    return(this_table)
        

def process_species(species, mode, basedir):
    path_els = [basedir, species, BUSCO_DATASET, "full_table.tsv"] 
    fulltablepath = "/".join(path_els)

    speciesTable = parse_table(fulltablepath, mode)
    filteredTable = speciesTable
    if mode == "transcriptome":
        filteredTable = filter_full_table(speciesTable)

    present_ids = list(filteredTable.keys())

    return (present_ids, filteredTable)

def process_all_species(species_list, basedir):
    retval = {} # {species: {table:{} , ids: {}}

    for species in species_list:
        this_ids, this_table = process_species(species, species_list[species], basedir)
        retval[species] = {"table": this_table, "ids": this_ids}
        del this_table
        del this_ids

    return retval

def process_all_fastas(species_list, busco_dict, basedir):
    retval = {} # {species: {table:{} , ids: {}}

    for species in species_list:
        this_fasta = parse_fasta(species, busco_dict, basedir)
        retval[species] = this_fasta
        del this_fasta

    for sp in retval:
        for busco in busco_dict:
            try:
                _ = retval[sp][busco]
            except KeyError:
                retval[sp][busco] = {"busco": busco, "gene": "", "seq": "", "prot": ""}

    return retval
    
    
def common_buscos(species_dict):

    buscos_list = [set(species_dict[x]["ids"]) for x in species_dict]
    common_buscos = reduce(lambda a,b: a.intersection(b), buscos_list)
    common_buscos = {k:{sp:species_dict[sp]["table"][k]["sequence"] for sp in species_dict} for k in common_buscos}

    return common_buscos


def fasta_from_dict(fasta_entry, seq_type="seq"):
    header = [">"+fasta_entry["species"], fasta_entry["busco"], fasta_entry["gene"]]
    header = "|".join(header)
    seq = fasta_entry[seq_type]
    retstr = header + "\n" + seq
    return retstr

def write_gene_fastas(fasta_dict, ouDir):
    '''
    This should iterate over the common genes found in all fastas and 
    write gene-wise fasta outputs with species names in the headers.

    e.g. a <gene_id>.fa file should have records like
    >Species_1
    SOMESEQUENCEEEEEE
    >Sepcies_2
    SOMEOTHERSEQUENCE
    '''

    common_final = reduce(lambda a,b: set(a).intersection(set(b)), [fasta_dict[sp].keys() for sp in fasta_dict])
    splist = list(fasta_dict.keys())
    depth = len(splist[0].split("/"))-1
    ##ouDir = "/".join(splist[0].split("/")[0:2]) + "/" + ouDir
    os.mkdir(ouDir)
    ##splist = [x.split("/")[depth] for x in splist]

    for busco in common_final:
        f = open(ouDir+ "/" + busco+'.fa', 'w')
        f2= open(ouDir+ "/" + busco+'.faa','w')
        for sp in splist:
            tmp = fasta_dict[sp][busco]
            tmp["species"] = sp.split("/")[depth]
            faseq = fasta_from_dict(tmp)
            faaseq= fasta_from_dict(tmp, "prot")
            f.write(faseq+'\n')
            f2.write(faaseq+'\n')
    return None

def write_gene_presence_report(species_dict, species_list, all_buscos, ouF):
    gene_presence = {"genome": [] , "transcriptome": []}
    denominators  = dict(Counter(list(species_list.values())))

    with open(ouF, "w") as f:
        ## First write a header
        f.write("Species")
        for busco in all_buscos:
            f.write("\t"+busco)
        f.write("\n")
        for species in species_dict:
            this_species = []
            f.write(species)
            for busco in all_buscos:
                status = int(busco in species_dict[species]["ids"])
                this_species.append(status)
                f.write("\t"+str(status))
            f.write("\n")
            mode = species_list[species]
            gene_presence[mode].append(this_species)

    ge_matrix = np.array(gene_presence["genome"]).transpose()
    tr_matrix = np.array(gene_presence["transcriptome"]).transpose()

    if ge_matrix.size > 0:
        ge_sums   = (ge_matrix.sum(axis=1) / denominators["genome"] * 100).round(decimals=2)
    else:
        ge_sums   = np.full((len(all_buscos),), 100)
    if tr_matrix.size > 0:
        tr_sums   = (tr_matrix.sum(axis=1) / denominators["transcriptome"] * 100).round(decimals=2)
    else:
        tr_sums   = np.full((len(all_buscos),), 100)
    genewise_report = {}
    for i in range(len(all_buscos)):
        tr_percent = tr_sums[i]
        ge_percent = ge_sums[i]
        gene = all_buscos[i]

        genewise_report[gene] = {"transcriptome": tr_percent,
                                 "genome": ge_percent}

    return (genewise_report, gene_presence)


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
        this_id = prot_alignment[sp_idx].id
        this_prot = str(prot_alignment[sp_idx,:].seq)
        this_prot_n = this_prot.replace(gap_char, "")

        this_nucl = codon_seqs[sp_idx].seq
        this_nucl_tra = this_nucl.upper().translate(table=codon_table)
        this_nucl = str(this_nucl)

        if this_nucl_tra != this_prot_n:
            ### Print an error/warning
            print(f"Error in alignment function at species {sp_idx}")
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

        this_record = SeqRecord(Seq(this_nucl_align), id=this_id.split("|")[0],
                                description="")

        nucl_records.append(this_record)

    return MultipleSeqAlignment(nucl_records)

def parse_alignment(gene):
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
    nucl_seqs = list(SeqIO.parse(nucl_unaligned_fname, "fasta"))
    codon_alignment = prot_to_codon(prot_alignment, nucl_seqs)

    if os.path.exists(gblocks_output_fname) or os.path.exists(guidance_output_fname) or os.path.exists(trimal_output_fname):
        try:
            clean_alignment = reduce(lambda a,b: a + b, [codon_alignment[:,x[0]:x[1]] for x in nucl_flanks])
        except TypeError:
            print(gene)
    else:
        clean_alignment = codon_alignment
    
    clean_alignment.sort()

    return clean_alignment

def process_all_alignments(alignment_output_dir):

    all_genes = [x.replace(".fa", "") for x in glob.glob(alignment_output_dir+"/*.fa")]
    cleaned_alignments = {gene.replace(alignment_output_dir+"/", ""):parse_alignment(gene) for gene in all_genes}

    return cleaned_alignments

def concat_alignments(alignment_dict):
    retval = reduce(lambda a,b: a + b, list(alignment_dict.values()))

    return retval

def combine_alignments(alignment_dict):
    iteron = [(k,Nexus.Nexus(format(alignment_dict[k], "nexus"))) for k in alignment_dict]
    combined_nexus = Nexus.combine(iteron)
    return(combined_nexus)
