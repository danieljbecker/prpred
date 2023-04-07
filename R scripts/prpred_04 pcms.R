## Plasmodium relictum prediction (prpred)
## 04_pcms
## danbeck@ou.edu
## last updated 4/6/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## libraries
library(ape)
library(caper)
library(phytools)
library(reshape2)
library(phylofactor)
library(ggtree)
library(treeio)

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

## trim to extant
data=data[data$Species.Status=="Extant",]
data$Species.Status=NULL

## clean
data=data[c("tip","TipLabel","Pr","lcites","Family3","Order3")]

## load BirdTree consensus phylogeny
setwd("/Users/danielbecker/Desktop/prpred/BirdTree phylo")
htree=readRDS("BirdNet 2K tree consensus.rds")

## fix labels
htree$tip.label=gsub("_"," ",htree$tip.label)

## save global tree
tree=htree

## combine data and phylogeny
data$label=data$tip
cdata=comparative.data(phy=tree,data=data,names.col=label,vcv=T,na.omit=F,warn.dropped=T)

## species
cdata$data$Species=cdata$data$tip

## phylogenetic signal in response for true and pseudo
## D of 0 = Brownian model, D of 1 = random (no phylogenetic signal)
set.seed(1)
mod=phylo.d(cdata,binvar=Pr,permut=1000)
mod ## D = 0.92, p < 0.001

## save label
cdata$data$label=cdata$data$tip

## get genus
cdata$data$Genus=sapply(strsplit(cdata$data$tip," "),function(x) x[1])

## taxonomy
cdata$data$taxonomy=paste(cdata$data$Order3,
                          cdata$data$Family3,
                          cdata$data$Genus,
                          cdata$data$Species,sep='; ')

## Holm rejection procedure
HolmProcedure <- function(pf,FWER=0.05){
  
  ## get split variable
  cs=names(coef(pf$models[[1]]))[-1]
  split=ifelse(length(cs)>2,cs[3],cs[1])
  
  ## obtain p values
  if (pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|t|)'])
  } else {
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|z|)'])
  }
  D <- length(pf$tree$tip.label)
  
  ## this is the line for Holm's sequentially rejective cutoff
  keepers <- pvals<=(FWER/(2*D-3 - 2*(0:(pf$nfactors-1))))
  
  
  if (!all(keepers)){
    nfactors <- min(which(!keepers))-1
  } else {
    nfactors <- pf$nfactors
  }
  return(nfactors)
}

## get species in a clade
cladeget=function(pf,factor){
  spp=pf$tree$tip.label[pf$groups[[factor]][[1]]]
  return(spp)
}

## summarize pf object 
pfsum=function(pf){
  
  ## get formula
  chars=as.character(pf$frmla.phylo)[-1]
  #chars=strsplit(as.character(pf$frmla.phylo)," ~ ")[[1]]
  
  ## response
  resp=chars[1]
  
  ## holm
  hp=HolmProcedure(pf)
  
  ## save model
  model=chars[2]
  
  ## set key
  setkey(pf$Data,'Species')
  
  ## make data
  dat=data.frame(pf$Data)
  
  ## make clade columns in data
  for(i in 1:hp){
    
    dat[,paste0(resp,'_pf',i)]=ifelse(dat$Species%in%cladeget(pf,i),'factor','other')
    
  }
  
  ## make data frame to store taxa name, response, mean, and other
  results=data.frame(matrix(ncol=6, nrow = hp))
  colnames(results)=c('factor','taxa','tips','node',"clade",'other')
  
  ## set taxonomy
  taxonomy=dat[c('Species','taxonomy')]
  taxonomy$taxonomy=as.character(taxonomy$taxonomy)
  
  ## loop
  for(i in 1:hp){
    
    ## get taxa
    tx=pf.taxa(pf,taxonomy,factor=i)$group1
    
    ## get tail
    tx=sapply(strsplit(tx,'; '),function(x) tail(x,1))
    
    ## combine
    tx=paste(tx,collapse=', ')
    
    # save
    results[i,'factor']=i
    results[i,'taxa']=tx
    
    ## get node
    tips=cladeget(pf,i)
    node=ggtree::MRCA(pf$tree,tips)
    results[i,'tips']=length(tips)
    results[i,'node']=ifelse(is.null(node) & length(tips)==1,'species',
                             ifelse(is.null(node) & length(tips)!=1,NA,node))
    
    ## get means
    ms=(tapply(dat[,resp],dat[,paste0(resp,'_pf',i)],mean))
    
    ## add in
    results[i,'clade']=ms['factor']
    results[i,'other']=ms['other']
    
  }
  
  ## return
  return(list(set=dat,results=results))
}

## binary, adjust for log1p citations
set.seed(1)
bpf=gpf(Data=cdata$data,tree=cdata$phy,
        frmla.phylo=Pr~phylo+lcites,
        family=binomial,algorithm='phylo',nfactors=4,min.group.size=20)

## summarize
bpf_results=pfsum(bpf)$results

## load in Pr edgelist
setwd("/Users/danielbecker/Desktop/prpred/MalAvi flat files")
pr_edge=read.csv("MalAvi edgelist_harmonized_Pr.csv")

## prune to extant
pr_edge=pr_edge[pr_edge$tip%in%data$tip,]

## simplify
pr_edge=pr_edge[c("tip","Lineage_Name")]

## trim to species in edgelist
htree=keep.tip(htree,unique(pr_edge$tip))

## load MalAvi phylogeny
setwd("/Users/danielbecker/Desktop/urbmigmalavi/MalAvi phylo")
ptree=read.tree("FastTree_output_tree.nhx")

## check correct match
spec=unique(pr_edge$Lineage_Name)
ptree=keep.tip(ptree,spec)
rm(spec)

## host breadth
table(pr_edge$Lineage_Name)

## make a H-P association matrix
hptab=table(pr_edge$tip,pr_edge$Lineage_Name)

## binary
hptab=ifelse(hptab>0,1,0)

# melt and visualize the matrix (computationally taxing)
hptab_melt=melt(hptab)
names(hptab_melt)=c("tip","Lineage_Name","value")
ggplot(hptab_melt,aes(tip,Lineage_Name))+
  geom_tile(aes(fill=value))+
  scale_fill_gradient(low="wheat1",high="steelblue")+
  theme_bw()+
  coord_flip()+
  theme(axis.text=element_blank(),
        axis.title=element_blank())+
  guides(fill=F)

## phylogenetic distances for host and parasite
htree_dist=cophenetic.phylo(htree)
ptree_dist=cophenetic.phylo(ptree)

## check
length(colnames(htree_dist))==length(unique(pr_edge$tip))
length(colnames(ptree_dist))==length(unique(pr_edge$Lineage_Name))

## paco
D=prepare_paco_data(H=htree_dist,P=ptree_dist,HP=hptab)
D=add_pcoord(D,correction="cailliez")
pac=PACo(D,nperm=999,seed=1,method="r0",symmetric=F) ## assumes column group tracks row group

## get interaction-specific cophylogenetic contributions
pac_links=paco_links(pac)

## pull out contributions
res=residuals_paco(pac_links$proc)
rdata=data.frame(res)
rdata$pairs=rownames(rdata)
x=strsplit(rdata$pairs,"-")
rdata$tip=sapply(x,function(x) x[1])
rdata$Lineage_Name=sapply(x,function(x) x[2])

## jackknife
jdata=data.frame(pac_links$jackknife)
names(jdata)=c("jackknife")
jdata$pairs=rownames(jdata)

## combine
rdata=merge(rdata,jdata,by="pairs",all.x=T)
rm(jdata)

## weight as rescaled res
wei=plotrix::rescale(res,c(0.5,0.05)) ## large res = little weight
#wei=((res^-2)/4)

## make the interaction matrix for cophyloplot format
imat=hptab_melt
imat=imat[which(imat$value>0),]
imat$value=NULL

## cophyloplot
par(oma=c(0,0,0,0),mar=c(0,0,0,0))
cophyloplot(htree,ptree,assoc=imat,show.tip.label=F,
            use.edge.length=F,lwd=wei,space=200,gap=5,length.line=0)
## lwd = wei