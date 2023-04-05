## Plasmodium relictum prediction (prpred)
## 01_MalAvi clean
## danbeck@ou.edu
## last updated 4/5/2023

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

## aggregate
hstab$id=with(hstab,paste(Lineage_Name,species,sep="_"))
hspr=aggregate(cbind(found,tested)~id,data=hstab,sum)

## bind in with hspr
hset=hstab[!duplicated(hstab$id),]
hset=hset[c("Lineage_Name","order","family","genus","species","id","Pr")]
hdata=merge(hset,hspr,by="id",all.x=T)

## clean
rm(hstab,hset,hspr)

## save
data=hdata
rm(hdata)

## flag hybrids
x=strsplit(data$species," ")
data$length=sapply(x,length)

## flag genus only
x=sapply(x,function(x) x[2])
data$sp=x

## remove hybrids and genus only
data=data[!data$length>2,]
data=data[!data$sp%in%c("sp","sp.","spp","spp."),]

## clean
data$length=NULL
data$sp=NULL

## match taxonomy
tax=taxonomy

## merge in with data
data=merge(data,tax,by="species",all.x=T)

## trim to host binary
bindata=data[!duplicated(data$species),]
bindata=bindata[c("species","id","order","family","genus","Pr","Jetz.species","match")]

## trim to Pr only
prdata=data[data$Pr==1,]

## clean
rm(lins,tax,data)

## export flat files
setwd("/Users/danielbecker/Desktop/prpred/MalAvi flat files")
write.csv(bindata,"binary Plasmodium relictum status MalAvi.csv")
write.csv(prdata,"Plasmodium relictum positivity MalAvi.csv")