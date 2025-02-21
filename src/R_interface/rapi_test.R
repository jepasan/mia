library(Rcpp)
library(mia)
library(miaSim)
library(ape)
library(picante)

source = "src/faith_R.cpp"
sourceCpp(source)

data(GlobalPatterns, package = "mia") # Matches flawlessly, no difference between picantes
data(esophagus, package = "mia") # Matches flawlessly, no difference between picantes
data(Tengeler2020, package = "mia") # Seems to work fine now, despite the tree being unrooted
tse <- GlobalPatterns

tse <- microbiomeDataSets::artificialgut() # All good (minor differences)
tse <- microbiomeDataSets::baboongut() # All good (minor differences)
tse <- microbiomeDataSets::SprockettTHData() # Matches picante with include.root=TRUE
tse <- microbiomeDataSets::GrieneisenTSData() # All good (minor differences)
tse <- microbiomeDataSets::qa10934() # Matches picante with include.root=TRUE
tse <- microbiomeDataSets::SilvermanAGutData() # All good (minor differences)

library(biomformat)
biom <- read_biom("/home/grunnar/downloads/con/finrisk.biom")
tree <- ape::read.tree("/home/grunnar/downloads/con/anonymized-finrisk-16S-BL_AGE.tre")
tse <- convertFromBIOM(biom)
rowTree(tse) <- tree
tse <- tse[1:100,1:100]

x <- addAlpha(tse, index="faith", fast_faith=TRUE)
faith <- colData(x)$faith
x2 <- addAlpha(tse, index="faith", fast_faith=FALSE)
faith2 <- colData(x2)$faith
x3 <- picante::pd(t(assay(tse)), rowTree(tse), include.root = TRUE)
faith3 <- as.vector(x3[[1]])
x4 <- picante::pd(t(assay(tse)), rowTree(tse), include.root = FALSE)
faith4 <- as.vector(x4[[1]])

faith - faith2
faith - faith3
faith - faith4

faith2 - faith3
faith3 - faith4

sum(abs(faith-faith2) > 0.00000001)
sum(abs(faith-faith3) > 0.00000001)
sum(abs(faith-faith4) > 0.00000001)


ape::write.tree(rowTree(tse), "//utuhome.utu.fi/jealpa/downloads/agut.tre")
write.biom()


#Checks
#This functions assumes that the tree representation is in cladewise order - Ensure this with ape's reorder.phylo() function
#Ensure that the tree is non-empty, etc
#Ensure that the tree doesn't get modified at any point
#Which assays are normally used for the calculations? - handled in the r functions

# The algorithm only calculates include.root = true values - These can differ significantly from the opposite case! 
# Possibly correctable in the c++ code?

y <- assay(tse)
t2 <- rowTree(tse)
f4 <- picante::pd(t(y), t2)[1]

faith - f4

rowTree(tse) <- ape::reorder.phylo(rowTree(tse), "cladewise")



picante::pd(t(ass), tree, include.root = FALSE)[1]


samples <- 10
obs <- 2000
tree <- ape::rtree(obs)
obsnames <- tree$tip.label
samplenames <- paste0("s", 1:samples)

#v <- rbinom(samples*obs, 1, 0.2) * rgeom(samples*obs, 0.05)
v <- rgeom(samples*obs, 0.05)

ass <- matrix(v, nrow=obs, ncol = samples)
colnames(ass) <- samplenames
rownames(ass) <- obsnames

testtree <- TreeSummarizedExperiment(assays=SimpleList(counts=ass))
rownames(testtree) <- obsnames
colnames(testtree) <- samplenames
rowTree(testtree) <- tree

z <- addAlpha(testtree, index="faith", fast_faith=TRUE)
faith_r <- colData(z)$faith
z2 <- addAlpha(testtree, index="faith", fast_faith=FALSE)
faith_r2 <- colData(z2)$faith
z3 <- picante::pd(t(assay(testtree)), rowTree(testtree), include.root = TRUE)
faith_r3 <- as.vector(z3[[1]])
z4 <- picante::pd(t(assay(testtree)), rowTree(testtree), include.root = FALSE)
faith_r4 <- as.vector(z4[[1]])


faith_r - faith_r2
faith_r - faith_r3
faith_r - faith_r4

faith_r2 - faith_r3
faith_r2 - faith_r4
faith_r3 - faith_r4

# So far, every randomly generated sample matches picante with include.root=TRUE
# in maybe 90% of cases the match the other methods, as well

# good samples for examining
#partiallyworks <- testtree # Matches picante with include.root=TRUE
#oldcodefails<- testtree # Seems like the old code fails with single-species samples in some cases?
            #Specifically when ape::drop.tip tries to reduce a two-edge tree to a single edge...
            #Otherwise matches picante with include.root=TRUE


edges <- c(1,2,3,4,5,7,8,10,11,12,14,19,20,21,23,24)
sum(tree$edge.length[edges])

# The current version of the code doesn't include a way to choose whether to include root or not
# Manually calculating the sum of edges for this tree produces a result that matches fastfaith!
# What is the old version doing differently?
# This seems to happen when the tree is simplified at the end of the process - collapse.singles eliminates the root node
# faith

.prune_tree <- function(treent, nodes){
    # Get those tips that can not be found from provided nodes
    remove_tips <- treent$tip.label[!treent$tip.label %in% nodes]
    # As long as there are tips to be dropped, run the loop
    while( length(remove_tips) > 0 ){
        # Drop tips that cannot be found. Drop only one layer at the time. Some
        # dataset might have taxa that are not in tip layer but they are in
        # higher rank. If we delete more than one layer at the time, we might
        # loose the node for those taxa. --> The result of pruning is a tree
        # whose all tips can be found provided nodes i.e., rows of TreeSE. Some
        # taxa might be higher rank meaning that all rows might not be in tips
        # even after pruning; these rows have still child-nodes that represent
        # other rows.
        # Suppress warning: drop all tips of the tree: returning NULL
        suppressWarnings(
            treent <- ape::drop.tip(
                treent, remove_tips,
                trim.internal = FALSE,
                collapse.singles = FALSE)
        )
        # If all tips were dropped, the result is NULL --> stop loop
        if( is.null(treent) ){
            warning("Pruning resulted to empty tree.", call. = FALSE)
            break
        }
        # Again, get those tips of updated tree that cannot be found from
        # provided nodes
        remove_tips <- treent$tip.label[!treent$tip.label %in% nodes]
    }
    # Simplify the tree structure. Remove nodes that have only single
    # descendant.
    if( !is.null(treent) && length(treent$tip.label) > 1 && ape::has.singles(treent) ){
        treent <- ape::collapse.singles(treent)
    }
    return(treent)
}
treent <- tree
nodes <- present[[7]]
.prune_tree(treent, nodes )


    # Gets vector where number represent nth sample
    ss <- seq_len(ncol(ass))
    
    # Repeats taxa as many times there are samples, i.e. get all the
    # taxa that are analyzed in each sample.
    taxa <- rep(rownames(ass), length(ss))
    
    # Gets those taxa that are present/absent in each sample.
    # Gets one big list that combines
    # taxa from all the samples.
    present_combined <- taxa[ ass[, ss] > 0 ]
    
    # Gets how many taxa there are in each sample. 
    # After that, determines indices of samples' first taxa with cumsum.
    split_present <- as.vector(cumsum(colSums(ass > 0)))
    
    # Determines which taxa belongs to which sample by first determining
    # the splitting points,
    # and after that giving every taxa number which tells their sample.
    split_present <- as.factor(cumsum((seq_along(present_combined)-1) %in%
                                          split_present))
    
    # Assigns taxa to right samples based on their number that they got from
    # previous step, and deletes unnecessary names.
    present <- unname(split(present_combined, split_present))
    
    # If there were samples without any taxa present/absent, the length of the
    # list is not the number of samples since these empty samples are missing.
    # Add empty samples as NULL.
    names(present) <- names(which(colSums2(ass) > 0))
    present[names(which(colSums2(ass) == 0))] <- list(NULL)
    present <- present[colnames(ass)]
    
    # Assign NA to all samples
    faiths <- rep(NA,length(ss))
    
    # If there are no taxa present, then faith is 0
    ind <- lengths(present) == 0
    faiths[ind] <- 0
    
    # If there are taxa present
    ind <- lengths(present) > 0
    # Loop through taxa that were found from each sample
    faiths_for_taxa_present <- lapply(present[ind], function(x){
        # Trim the tree
        temp <- .prune_tree(tree, x)
        # Sum up all the lengths of edges
        temp <- sum(temp$edge.length)
        return(temp)
    })
    faiths_for_taxa_present <- unlist(faiths_for_taxa_present)
    faiths[ind] <- faiths_for_taxa_present
    