from importlib import reload
import find_buscos_fragmented as find_buscos
import sys
import os
import yaml
reload(find_buscos)

'''
single folder run for bostest
'''
BUSCO_DATASET = "run_lepidoptera_odb10"
run_folder = sys.argv[1]
##prev_run = clade_list[0]

# fallback = "papilionoidea"
# fallback_keys = set(find_buscos.detect_busco_mode(find_buscos.detect_species_subfolders(fallback+"/busco_results")).keys())
# fallback_keys = set([x.split("/")[-1] for x in fallback_keys])

splist = {}
##species_list = find_buscos.detect_species_subfolders(run_folder+"/busco_results")
splist = find_buscos.detect_busco_mode(run_folder+"/busco_results")
##splist.update(this_splist)


print("Detected following species with the specified modes:")
print(yaml.dump(splist, default_flow_style=False))
print("====================================================")
#transthresh = 80 ; genomethresh = 90
transthresh = int(sys.argv[2]) ; genomethresh = int(sys.argv[3])
all_species = find_buscos.process_all_species(splist, run_folder+"/busco_results")
allBuscos   = find_buscos.parse_table(run_folder+"/busco_results/" + list(splist.keys())[0] +"/"+ BUSCO_DATASET+ "/full_table.tsv", "ids")
print("Calculating gene presences across all species...")
## genewise, gene_report = find_buscos.write_gene_presence_report(all_species, splist, allBuscos, "gene_presence.tsv") 
genewise, gene_report = find_buscos.write_gene_presence_report(all_species, splist, allBuscos, run_folder+"/gene_presence.tsv") 
passedGenes = [x for x in genewise if (genewise[x]["transcriptome"] >= transthresh and genewise[x]["genome"] >= genomethresh)]
print("Done.")
print(f"{len(passedGenes)} Genes are present in at least %{genomethresh} of the genomes and %{transthresh} of the transcriptomes")

passed_genes = {k:{sp:'' for sp in all_species} for k in passedGenes}
for k in passedGenes:
    for sp in all_species:
        try:
            passed_genes[k][sp] = all_species[sp]["table"][k]["rec_id"]
        except KeyError:
            passed_genes[k][sp] = ''

# processed = [x for x in os.listdir(prev_run+"/mafft_output_new_clades") if x.endswith(".faa")]
# [processed.append(x) for x in os.listdir(prev_run+"/mafft_output_never_processed") if x.endswith(".faa")]
# ## By default, "mafft_output_new_clades" folder in the "prev_run" folder
# ## should include alignments from the previous run when that folder was
# ## holding output and information as the added_clade of that run

# processed = [x.replace("_aligned.faa", "") for x in processed]
# processed = set(processed)
# pg = set(list(passed_genes.keys()))
# remainder_genes = pg.difference(processed) ## These were not processed during the ref clade's prev run

# pg_processed = {}         ## gene information for the new clades
# pg_never_processed = {}     ## gene information of the new genes for the ref clade

# for busco in passed_genes:
#     if busco in remainder_genes:
#         this_never_processed = {k:passed_genes[busco][k] for k in passed_genes[busco]}
#         pg_never_processed[busco] = this_never_processed
#     else:
#         this_processed = {k:passed_genes[busco][k] for k in passed_genes[busco] if not k.startswith(prev_run)}
#         pg_processed[busco] = this_processed

# splist_ref = {k:splist[k] for k in splist if k.startswith(prev_run)} #and k.split("/")[2] in papilionoidea_remainder}
# splist_f = {k:splist[k] for k in splist if not k.startswith(prev_run)}
# splist_no_fallback = {k:splist[k] for k in splist if k.split("/")[-1] not in fallback_keys }

print("Parsing fastas to search for selected genes...")
fastas = find_buscos.process_all_fastas(splist, passed_genes , run_folder+"/busco_results")
print("done.")

print("Writing gene-wise codon and protein fastas...")
find_buscos.write_gene_fastas(fastas, run_folder+"/genewise_fastas")
print("done.")

