library(Rcpp)
library(mia)
library(biomformat)
library(ape)
library(rhdf5)

equals <- function(x, y, msg){
    if (x!=y)
        stop(msg)


}
aboutEquals <- function(x, y, msg){
	if((x-y)>0.005)
        stop(msg)
}

source = "R/unifrac_cpp/su_R_s.cpp"
sourceCpp(source)
table = "test.biom"
tree = "test.tre"
nthreads = 1

b <- read_hdf5_biom("R/unifrac_cpp/R_interface/test.biom")
outfile <- tempfile()
write_biom(b, outfile)
bb <- read_biom(outfile)


fname <- "R/unifrac_cpp/R_interface/test.tre"
tree <- ape::rtree(500)
ape::write.tree(tree, file = fname, append = FALSE)

newick <- readChar(fname, file.info(fname)$size)
tree <- ape::read.tree(fname)

identical(treetest(newick), treetest2(tree))
z <- treetest(newick)

tree2 <- ape::reorder.phylo(tree, "postorder")
tree3 <- ape::read.tree("R/unifrac_cpp/R_interface/test2.tre")

treese <- makeTreeSEFromBiom(bb, treefilename=tree)
treese2 <- changeTree(treese, rowTree = tree)


data(GlobalPatterns, package = "mia")
data(esophagus, package = "mia")
tse <- esophagus


#faith = faith_pd(table, tree)
faith_pd_new(tse)



assays(tse)[[1]]
colSums(assays(tse)[[1]])

exp = c(4, 5, 6, 3, 2, 5)

equals(faith["n_samples"][[1]], 6, "n_samples != 6")
for ( i in 1:6){
	aboutEquals(faith["faith_pd"][[1]][i], exp[i], "Output not as expected")
}

print('Success.')

print('All tests pass')


#Checks
#This functions assumes that the tree representation is in cladewise order - Ensure this with ape's reorder.phylo() function
#Ensure that the tree is non-empty, etc
#Ensure that the tree doesn't get modified at any point
#Which assays are normally used for the calculations?