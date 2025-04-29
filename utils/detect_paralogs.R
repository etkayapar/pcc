library(ape)
args = commandArgs(trailingOnly=TRUE)
## args = c(0.08, 0.55,
##          "../output/after_trimal/outlier_detection/saved_genes_removed_taxa/output.treefile",
##          "../output/after_trimal/outlier_detection/all_genes_names.txt",
##          "./tmp_output")
print(args)
threshold = args[1]
tax_threshold = args[2]
trees_file = args[3] ## Should point to gene_trees file output by TreeShrink
gene_names_file = args[4]
output_dir=args[5]

## Functions ------
get_internal_relbrlens = function(tree){
    relbrlens = tree$edge.length / sum(tree$edge.length)
    internalbr_ids = which(tree$edge[,2] > Ntip(tree))
    internal_relbrlens = relbrlens[internalbr_ids]
    return(internal_relbrlens)
}

prop.max.f = function(tree){
    internal_relbrlens = get_internal_relbrlens(tree)
    return(max(internal_relbrlens))
}

split_trees = function(tree, retain_threshold=0.9, prop_max_threshold=0.1){
    require(phytools)
    ##p_m_t = prop_max_threshold
    prop_max = prop.max.f(tree)
    relbrlens = tree$edge.length / sum(tree$edge.length)
    ## Check if there is a longer than "threshold" prop_max branch in the tree
    ## after splitting. If there is, re-run the function until it goes away
    if(prop_max <= prop_max_threshold){
        print(prop_max)
        return(tree)
    }else{
        cat("one round of splitting...")
    }
    ntips = Ntip(tree)
    longest_branch = which(relbrlens == prop_max)
    rootwards = tree$edge[longest_branch,1]
    tipwards = tree$edge[longest_branch, 2]

    splits = splitTree(tree, list(node=tipwards, bp=1))
    rootwards_split = splits[[1]]
    rootwards_split = drop.tip(rootwards_split, tip="NA")
    tipwards_split = splits[[2]]
    rootwards_ptips = Ntip(rootwards_split) / ntips
    tipwards_ptips = Ntip(tipwards_split) / ntips

    if(tipwards_ptips > rootwards_ptips && tipwards_ptips >= retain_threshold){
        cat("splitting tipwards...")
        return(split_trees(tipwards_split, prop_max_threshold = prop_max_threshold))
    }else if(rootwards_ptips > tipwards_ptips && rootwards_ptips >=retain_threshold){
        cat("splitting rootwards...")
        return(split_trees(rootwards_split, prop_max_threshold = prop_max_threshold))
    }else{
        return(NULL)
    }
}

### Reading in the tree and the data -----
gene_trees = read.tree(trees_file)
gene_names = scan(gene_names_file, what=character())
#setwd(output_dir)

#max_n_tax = max(sapply(gene_trees, function(x){length(x$tip.label)}))
max_n_tax = length(unique(unlist(sapply(gene_trees, function(x){return(x$tip.label)}))))

### Calculating the relative length of the longest branch ----
prop_max = sapply(gene_trees, prop.max.f)
names(prop_max) = gene_names
names(gene_trees) = gene_names

prop_tax = sapply(gene_trees, function(x){length(x$tip.label)/max_n_tax})
names(prop_tax) = gene_names
### Filtering and processing trees according to defined thresholds ----

long_branch_enough_taxa_genes = names(prop_max[prop_max > threshold])

#### Running tree splitting function on suitable genes -----

saved_trees = lapply(long_branch_enough_taxa_genes, function(gene){
    this_tree = gene_trees[[gene]]
    this_tree_saved = split_trees(this_tree, retain_threshold = tax_threshold, prop_max_threshold = threshold)
    if(!is.null(this_tree_saved)){
        message(paste(gene, "saved"))
        if(Ntip(this_tree_saved)/Ntip(this_tree) < tax_threshold){
            message(paste(gene, "saved but has not enough taxa"))
            this_tree_saved = NULL
        }
    }
    return(this_tree_saved)
})
## message(paste(length(saved_trees), "trees were saved"))
## message(paste(sum(sapply(saved_trees, is.null)), "trees were lost among them"))
names(saved_trees) = long_branch_enough_taxa_genes
pdf(paste(output_dir,"saved_genes.pdf", sep="/"), 40,20)
if(length(saved_trees) !=0) {
    print(str(saved_trees))
    print(summary(saved_trees))
    par(mfrow=c(1,2))
    sapply(names(saved_trees), function(gene){
        orig_tree = gene_trees[[gene]]
        saved_tree = saved_trees[[gene]]
        plot(orig_tree, cex=0.5, main=gene)
        if(is.null(saved_tree)){
            plot.new()
        }else{
            plot(saved_tree, cex=0.5, main=paste(gene,"_saved", sep=""))
        }
        return(NULL)
    })
}
dev.off()


if(length(saved_trees) !=0) {
    saved_genes = names(saved_trees)[sapply(saved_trees, function(x){!is.null(x)})]
}

#### Write out kept taxa for saved genes ------
    out.dir = output_dir
    ## if (!file.exists(out.dir)){
    ##     dir.create(out.dir)
    ## }

if(length(saved_trees) !=0) {
    sapply(saved_genes, function(gene){
        nm = paste(gene, "_removed_taxa_r.txt", sep="")
        nm = file.path(out.dir, nm)
        kept_taxa = saved_trees[[gene]]$tip.label
        removed_taxa = setdiff(gene_trees[[gene]]$tip.label, kept_taxa)
        write(removed_taxa, file=nm, sep='\n')
    })
}
