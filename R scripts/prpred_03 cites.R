## Plasmodium relictum prediction (prpred)
## 03_cites
## danbeck@ou.edu
## last updated 4/1/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## libraries
library(easyPubMed)

## load in MalAvi with traits
setwd("/Users/danielbecker/Desktop/prpred")
data=read.csv("MalAvi with host traits.csv")

## collect citations per host species
cites=c()
for(i in 1:length(data$tip)) {
  
  ## split into genus and species
  x=strsplit(data$tip[i]," ")[[1]]
  
  ## paste
  x=paste(x[1],"AND",x[2])
  
  ## get citations
  counts=as.numeric(as.character(get_pubmed_ids(x)$Count))
  cites[i]=counts
  print(paste(i,"/",nrow(data)))
}

## compile all citations
cdata=data.frame(tip=data$tip,
                 cites=cites)



