## Plasmodium relictum prediction (prpred)
## 04_pcms
## danbeck@ou.edu
## last updated 4/5/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## libraries
library(ape)
library(caper)
library(phytools)

## load in MalAvi with traits
setwd("/Users/danielbecker/Desktop/prpred")
data=read.csv("MalAvi with host traits.csv")

## combine with cites
cdata=read.csv("BirdTree citations.csv")
data=merge(data,cdata,by="tip")
rm(cdata)

## clean
data$X.x=NULL
data$X.y=NULL

## fix citations
data$lcites=log1p(data$cites)
data$cites=NULL

## tabulate labels
table(data$Pr)
round(prop.table(table(data$Pr)),2)

## check extinct/invalid
table(data$Species.Status,data$Pr)

## trim to extant
data=data[data$Species.Status=="Extant",]
data$Species.Status=NULL

## clean
data=data[c("tip","TipLabel","Pr","lcites","Family3","Order3")]

## load phylogenies
