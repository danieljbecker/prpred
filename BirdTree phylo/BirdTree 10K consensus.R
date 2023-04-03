## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## libraries
library(ape)
library(caper)

### Combine Hackett Stage 2 trees (10K)
setwd("~/Desktop/BirdTree Hackett Stage 2/")
tree1=read.tree("mnt/data/projects/birdphylo/Tree_sets/Stage2_full_data/CombinedTrees/AllBirdsHackett1.tre")
tree2=read.tree("mnt 2/data/projects/birdphylo/Tree_sets/Stage2_full_data/CombinedTrees/BirdzillaHackett2.tre")
tree3=read.tree("mnt 3/data/projects/birdphylo/Tree_sets/Stage2_full_data/CombinedTrees/BirdzillaHackett3.tre")
tree4=read.tree("mnt 4/data/projects/birdphylo/Tree_sets/Stage2_full_data/CombinedTrees/BirdzillaHackett4.tre")
tree5=read.tree("mnt 5/data/projects/birdphylo/Tree_sets/Stage2_full_data/CombinedTrees/BirdzillaHackett5.tre")
tree6=read.tree("mnt 6/data/projects/birdphylo/Tree_sets/Stage2_full_data/CombinedTrees/BirdzillaHackett6.tre")
tree7=read.tree("mnt 7/data/projects/birdphylo/Tree_sets/Stage2_full_data/CombinedTrees/BirdzillaHackett7.tre")
tree8=read.tree("mnt 8/data/projects/birdphylo/Tree_sets/Stage2_full_data/CombinedTrees/BirdzillaHackett8.tre")
tree9=read.tree("mnt 9/data/projects/birdphylo/Tree_sets/Stage2_full_data/CombinedTrees/BirdzillaHackett9.tre")
tree10=read.tree("mnt 10/data/projects/birdphylo/Tree_sets/Stage2_full_data/CombinedTrees/BirdzillaHackett10.tre")

## combien into big multiPhylo
trees=c(tree1,tree2,tree3,tree4,tree5,tree6,tree7,tree8,tree9,tree10)

## randomly sample 2K trees
set.seed(1)
tsample=sample(trees,2000)

## make consensus
ctree=consensus(tsample)

## resolve
tree=multi2di(ctree)
tree=compute.brlen(tree)
tree=makeLabel(tree)

### Export consensus tree as is
setwd("~/OneDrive - University of Oklahoma/Becker Lab/Datasets/BirdTree")
saveRDS(tree,"BirdNet 2K tree consensus.rds")