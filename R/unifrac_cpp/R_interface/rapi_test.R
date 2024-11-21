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

source = "R/unifrac_cpp/su_R.cpp"
sourceCpp(source)
table = "test.biom"
tree = "test.tre"
nthreads = 1

b <- read_hdf5_biom("R/unifrac_cpp/R_interface/test.biom")
outfile <- tempfile()
write_biom(b, outfile)
bb <- read_biom(outfile)

fname <- "R/unifrac_cpp/R_interface/test.tre"
newick <- readChar(fname, file.info(fname)$size)
z <- treetest(newick)

tree <- ape::read.tree("R/unifrac_cpp/R_interface/test.tre")
treese <- makeTreeSEFromBiom(bb, treefilename=tree)

treese2 <- changeTree(treese, rowTree = tree)

test <- tempfile()
rowTree(treese2)

print('Testing Faith PD..')

#faith = faith_pd(table, tree)
faith <- faith_pd_new(treese2, newick)

exp = c(4, 5, 6, 3, 2, 5)

equals(faith["n_samples"][[1]], 6, "n_samples != 6")
for ( i in 1:6){
	aboutEquals(faith["faith_pd"][[1]][i], exp[i], "Output not as expected")
}

print('Success.')

print('All tests pass')

