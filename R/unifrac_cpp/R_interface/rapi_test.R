library(Rcpp)
library(mia)
library(miaSim)
library(ape)
library(picante)

source = "R/unifrac_cpp/su_R.cpp"
sourceCpp(source)

data(GlobalPatterns, package = "mia")
data(esophagus, package = "mia")
data(HintikkaXOData, package = "mia")
data(Tengeler2020, package = "mia") # This dataset produces divergent values for some reason - Presumably something to do with the tree being unrooted
tse <- GlobalPatterns
rowTree(tse) <- ape::reorder.phylo(rowTree(tse), "cladewise")

ts1 <- rowTree(tse)

fname <- "R/unifrac_cpp/R_interface/tree.tre"
ape::write.tree(ts1, fname)

newick <- readChar(fname, file.info(fname)$size)

y <- rowTree_to_bp(ts1)
x <- newick_to_bp(newick)

faith <- faith_cpp(tse)
x <- estimateDiversity(tse, index="faith")
faith2 <- colData(x)$faith

faith - faith2

#Checks
#This functions assumes that the tree representation is in cladewise order - Ensure this with ape's reorder.phylo() function
#Ensure that the tree is non-empty, etc
#Ensure that the tree doesn't get modified at any point
#Which assays are normally used for the calculations?





tse2 <- estimateDiversity(tse_hubbell, index = "faith")
colData(tse2)$faith
faith_pd(tse_hubbell)


t <- ape::rtree(12706, rooted = F, tip.label = rownames(tse2))
tse2 <- tse[[1]]
rowTree(tse2) <- t
tse2 <- estimateDiversity(tse2)
faith2 <- colData(tse2)$faith
faith <- faith_pd(tse2, T)
faith - faith2

data(phylocom, package = "picante")
x <- phylocom$sample
y <- assay(tse)
t1 <- phylocom$phylo
t2 <- rowTree(tse)
picante::pd(x, t1)
picante::pd(t(y), t2, include.root=F)


rowTree(tse) <- ape::reorder.phylo(rowTree(tse), "cladewise")
