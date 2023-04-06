## Plasmodium relictum prediction (prpred)
## 01_MalAvi clean
## danbeck@ou.edu
## last updated 4/6/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## libraries
library(malaviR)

## obtain most recent lineage summary
lineage=extract_table("Grand Lineage Summary")

## extract all lineages for Plasmodium relictum
lset=lineage[lineage$species=="Plasmodium relictum",]

## get lineages
lins=lset$Lineage_Name

## clean
rm(lineage,lset)

## get host and sites table for full associations
hstab=extract_table("Hosts and Sites Table")

## make Pr
hstab$Pr=ifelse(hstab$Lineage_Name%in%lins,1,0)

## just edgelist
#edge=hstab[hstab$Pr==1,]
edge=hstab

## unique
edge$id=with(edge,paste(Lineage_Name,species,sep="_"))
edge$record=1

## aggregate
esum1=data.frame(aggregate(record~id,data=edge,sum))
esum2=data.frame(aggregate(Pr~id,data=edge,sum))
esum=merge(esum1,esum2,by="id")
rm(esum1,esum2)

## fix
esum$Pr=ifelse(esum$Pr>0,1,0)

## bind in with edge
edge=edge[!duplicated(edge$id),]
edge=edge[c("Lineage_Name","parasiteGenus","order","family","genus","species","id")]
edge=merge(edge,esum,by="id",all.x=T)

## clean
rm(hstab,lins,esum)

## flag hybrids
x=strsplit(edge$species," ")
edge$length=sapply(x,length)

## flag genus only
x=sapply(x,function(x) x[2])
edge$sp=x

## remove hybrids and genus only
edge=edge[!edge$length>2,]
edge=edge[!edge$sp%in%c("sp","sp.","spp","spp."),]

## clean
edge$length=NULL
edge$sp=NULL

## match taxonomy
tax=taxonomy

## merge in with data
edge=merge(edge,tax,by="species",all.x=T)

## clean
rm(x,tax)

## export flat file
setwd("/Users/danielbecker/Desktop/prpred/MalAvi flat files")
write.csv(edge,"MalAvi edgelist_cleaned.csv")