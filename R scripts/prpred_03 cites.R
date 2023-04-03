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

## collect citations per host species, using latin AND common
cites=c()
for(i in 1:length(data$tip)) {
  
  ## split into genus and species
  x=strsplit(data$tip[i]," ")[[1]]
  
  ## paste
  x=paste(x,collapse=" AND ")
  
  ## repeat for common name
  y=strsplit(data$English[i]," ")[[1]]
  y=paste(y,collapse=" AND ")
  
  ## add parentheses
  x=paste("(",x,")",sep="")
  y=paste("(",y,")",sep="")
  
  ## combine
  xy=paste(x," AND ",y,sep="")
  
  ## get citations
  counts=as.numeric(as.character(get_pubmed_ids(xy)$Count))
  cites[i]=counts
  print(paste(data$tip[i],", cites = ",counts,", ",i,"/",nrow(data),sep=""))
}

## compile all citations
cdata=data.frame(tip=data$tip,
                 cites=cites)

## export
setwd("/Users/danielbecker/Desktop/prpred")
write.csv(cdata,"BirdTree citations.csv")