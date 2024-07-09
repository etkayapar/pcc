library(ape)
args = commandArgs(trailingOnly=TRUE)
print(args)
threshold = args[1]
tax_threshold = args[2]
trees_file = args[3]
gene_names_file = args[4]
output_dir=args[5]

## Functions ------
prop.max.f = function(tree){
    return(max(tree$edge.length) / sum(tree$edge.length))
}
split_trees = function(tree, retain_threshold=0.9, prop_max_threshold=0.1){
    require(phytools)
    ##p_m_t = prop_max_threshold
    prop_max = prop.max.f(tree)
    ## Check if there is a longer than 0.1 prop_max branch in the tree
    ## after splitting. If there is, re-run the function until it goes away
    if(prop_max <= prop_max_threshold){
        print(prop_max)
        return(tree)
    }else{
        cat("one round of splitting...")
    }
    ntips = Ntip(tree)
    longest_branch = which.max(tree$edge.length)
    rootwards = tree$edge[longest_branch,1]
    tipwards = tree$edge[longest_branch, 2]
    
    splits = splitTree(tree, list(node=tipwards, bp=1))
    rootwards_split = splits[[1]]
    rootwards_split = drop.tip(rootwards_split, tip="NA")
    tipwards_split = splits[[2]]
    rootwards_ptips = Ntip(rootwards_split) / ntips
    tipwards_ptips = Ntip(tipwards_split) / ntips
    
    if(tipwards_ptips >= retain_threshold){
        cat("splitting tipwards...")
        return(split_trees(tipwards_split, prop_max_threshold = prop_max_threshold))
    }else if(rootwards_ptips >=retain_threshold){
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

max_n_tax = max(sapply(gene_trees, function(x){length(x$tip.label)}))

### Calculating the relative length of the longest branch ----
brlens = lapply(gene_trees,function(X){X$edge.len})
prop_max = sapply(brlens, function(X){max(X)/sum(X)})
names(prop_max) = gene_names
names(gene_trees) = gene_names

prop_tax = sapply(gene_trees, function(x){length(x$tip.label)/max_n_tax})
names(prop_tax) = gene_names
### Filtering and processing trees according to defined thresholds ----

long_branch_genes = prop_max[prop_max > threshold]
long_branch_enough_taxa_genes = names(prop_max[prop_max > threshold & prop_tax >= tax_threshold])

#### Running tree splitting function on suitable genes -----

saved_trees = sapply(long_branch_enough_taxa_genes, function(gene){
    this_tree = gene_trees[[gene]]
    this_tree_saved = split_trees(this_tree, retain_threshold = 0.95, prop_max_threshold = 0.1)
    if(!is.null(this_tree_saved)){
        if(Ntip(this_tree_saved)/Ntip(this_tree) < 0.95){
            this_tree_saved = NULL
        }
    }
    return(this_tree_saved)
})
pdf(paste(output_dir,"saved_genes.pdf", sep="/"), 40,20)
if(length(saved_trees) !=0) {
    sum(sapply(saved_trees, is.null))
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

outlier_genes = names(prop_max[prop_max > threshold | prop_tax < tax_threshold])

if(length(saved_trees) !=0) {
    outlier_genes = setdiff(outlier_genes, saved_genes)
}

pdf(paste(output_dir,"outlier_genes.pdf", sep="/"),30,30)
par(mfrow=c(5,5))
sapply(outlier_genes, function(X){
    plot(gene_trees[[X]], main=X)
})
dev.off()
write(outlier_genes, file = paste(output_dir,"outlier_genes.txt", sep="/"),sep='\n')

#### Write out kept taxa for saved genes ------
    out.dir = paste(output_dir,"saved_genes_kept_taxa", sep="/")
    if (!file.exists(out.dir)){
        dir.create(out.dir)
    }

if(length(saved_trees) !=0) {
    sapply(saved_genes, function(gene){
        nm = paste(gene, "_kept_taxa.txt", sep="")
        nm = file.path(output_dir,out.dir, nm)
        taxa = saved_trees[[gene]]$tip.label
        write(taxa, file=nm, sep='\n')
    })
}
