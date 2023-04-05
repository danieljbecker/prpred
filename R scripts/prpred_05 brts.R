## Plasmodium relictum prediction (prpred)
## 04_brts
## danbeck@ou.edu
## last updated 4/2/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## libraries
library(gbm)
library(fastDummies)
library(rsample)
library(ROCR)
library(sciplot)
library(ggplot2)
library(InformationValue)
library(PresenceAbsence)

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

## get families and orders of (screened or) positive birds
pos=data[data$Pr==1 | data$zero=="true zero",]
pos=data[data$Pr==1,]
porder=unique(pos$Order3) ## AVONET
pfam=unique(pos$Family3) ## AVONET

## positive coverage
length(porder)/length(unique(data$Order3))
length(pfam)/length(unique(data$Family3))

## tabulate when trimming by order
round(prop.table(table(data[data$Order3%in%porder,"Pr"])),2)
round(prop.table(table(data[data$Family3%in%pfam,"Pr"])),2)

## trim
#data=data[data$Order3%in%porder,]
data=data[data$Family3%in%pfam,]

## tabulate labels
table(data$Pr)
prop.table(table(data$Pr))

## clean traits
data$X=NULL
data$WJSpecID=NULL
data$Scientific=NULL
data$PatchClade=NULL
data$FileName=NULL
data$Taxo=NULL
data$BLFamilyEnglish=NULL
data$FamSequID=NULL
data$IOCOrder=NULL
data$PassNonPass=NULL
data$OscSubOsc=NULL
data$species=NULL
data$id=NULL
data$order=NULL
data$family=NULL
data$genus=NULL
data$Jetz.species=NULL
data$match=NULL
data$Total.individuals=NULL
data$Female=NULL
data$Male=NULL
data$Unknown=NULL
data$Complete.measures=NULL
data$Mass.Source=NULL
data$Reference.species=NULL
data$Mass.Refs.Other=NULL
data$Inference=NULL
data$Traits.inferred=NULL

## check extinct/invalid
table(data$Species.Status,data$Pr)

## trim to extant
data=data[data$Species.Status=="Extant",]
data$Species.Status=NULL

## clean IUCN
data$IUCN=ifelse(data$X2010.IUCN.Red.List.category=="",NA,data$X2010.IUCN.Red.List.category)
data$X2010.IUCN.Red.List.category=NULL

## tabulate class and number of levels
vars=melt(lapply(data, class))
names(vars)=c("class","var")
uset=data.frame(apply(data,2,function(x) length(unique(x))))
names(uset)=c("levels")
uset$var=rownames(uset)
vars=merge(vars,uset,by="var")
rm(uset)

## if categorical & levels > 5, binary
rm(vars)

## save raw
raw=data

## make binary columns for IUCN
dums=dummy_cols(raw["IUCN"])

## unique
dums=dums[!duplicated(dums$IUCN),]

## ensure all factor
for(i in 1:ncol(dums)){
  
  ## column as factor
  dums[,i]=factor(dums[,i])
  
}

## merge
data=merge(data,dums,by="IUCN",all.x=T)
rm(dums)
data$IUCN=NULL
data$IUCN_NA=NULL

## make binary columns for Habitat
dums=dummy_cols(raw["Habitat"])

## unique
dums=dums[!duplicated(dums$Habitat),]

## ensure all factor
for(i in 1:ncol(dums)){
  
  ## column as factor
  dums[,i]=factor(dums[,i])
  
}

## merge
data=merge(data,dums,by="Habitat",all.x=T)
rm(dums)
data$Habitat=NULL
data$Habitat_NA=NULL

## make binary for Trophic.Niche
dums=dummy_cols(raw["Trophic.Niche"])

## unique
dums=dums[!duplicated(dums$Trophic.Niche),]

## ensure all factor
for(i in 1:ncol(dums)){
  
  ## column as factor
  dums[,i]=factor(dums[,i])
  
}

## merge
data=merge(data,dums,by="Trophic.Niche",all.x=T)
rm(dums)
data$Trophic.Niche=NULL
data$Trophic.Niche_NA=NULL

## make binary columns for family
dums=dummy_cols(raw["Family3"])

## unique
dums=dums[!duplicated(dums$Family3),]

## ensure all factor
for(i in 1:ncol(dums)){
  
  ## column as factor
  dums[,i]=factor(dums[,i])
  
}

## merge
data=merge(data,dums,by="Family3",all.x=T)
rm(dums)
data$Family3=NULL

## repeat for order
dums=dummy_cols(raw["Order3"])

## unique
dums=dums[!duplicated(dums$Order3),]

## ensure all factor
for(i in 1:ncol(dums)){
  
  ## column as factor
  dums[,i]=factor(dums[,i])
  
}

## merge
data=merge(data,dums,by="Order3",all.x=T)
rm(dums)
data$Order3=NULL

## assess missing values
mval=data.frame(apply(data,2,function(x) length(x[!is.na(x)])/nrow(raw)))

## get names
mval$variables=rownames(mval)
names(mval)=c("comp","column")

## if have at least 50% values, keep
mval$comp=round(mval$comp,3)
mval$keep=ifelse(mval$comp>=0.50,"keep","cut")
mval=mval[order(mval$keep),]
keeps=mval[-which(mval$keep=="cut"),]$column

## order
mval=mval[order(mval$comp),]

## drop if not well represented
data=data[keeps]
rm(mval,keeps,porder,i)

## set
set=data
set$TipLabel=NULL
set$Hackett_FineClades=NULL
set$Hackett_CoarseClades=NULL
set$English=NULL
set$BLFamilyLatin=NULL
set$Species3=NULL
set$zero=NULL

## tally
ncol(set)-2 ## tip and Pr

## factors
set$Migration=factor(set$Migration)
set$Trophic.Level=factor(set$Trophic.Level)
set$Primary.Lifestyle=factor(set$Primary.Lifestyle)
set$urban=factor(set$urban)
set$humanDisturbed=factor(set$humanDisturbed)

## extract continuous vars
vars=names(which(unlist(lapply(set,class))!="factor"))
vars=vars[!vars%in%c("tip","Pr")]

## correlate
cmat=round(cor(set[vars],use="complete.obs"),2)
diag(cmat)=NA
cset=as.data.frame(as.table(cmat))
cset[cset$Freq>0.9,]

## remove min/max latitude
set$Min.Latitude=NULL
set$Max.Latitude=NULL

## remove secondary
set$Secondary1=NULL

## remove beak length nares
set$Beak.Length_Nares=NULL

## clean
rm(vars,pos,cset,cmat)

## tally
ncol(set)-2 ## tip and Pr

# ## hyperparameter grid
# hgrid=expand.grid(n.trees=6000,
#                   interaction.depth=c(2,3,4),
#                   shrinkage=c(0.01,0.001),
#                   n.minobsinnode=c(4,5),
#                   seed=seq(1,5,by=1),
#                   prop=c(0.8,0.9))
# 
# ## trees, depth, shrink, min, prop
# hgrid$id=with(hgrid,paste(n.trees,interaction.depth,shrinkage,n.minobsinnode))
# 
# ## sort by id then seed
# hgrid=hgrid[order(hgrid$id,hgrid$seed),]
# 
# ## now add rows
# hgrid$row=1:nrow(hgrid)
# 
# ## factor id
# hgrid$id2=factor(as.numeric(factor(hgrid$id)))
# 
# ## function to assess each hyperpar combination
# hfit=function(row){
#   
#   ## make new data
#   ndata=set
#   ndata$tip=NULL
#   
#   ## correct response
#   ndata$response=ndata$Pr
#   ndata$Pr=NULL
#   
#   ## use rsample to split
#   set.seed(hgrid$seed[row])
#   split=initial_split(ndata,prop=hgrid$prop[row],strata="response")
#   
#   ## test and train
#   dataTrain=training(split)
#   dataTest=testing(split)
#   
#   ## yTest and yTrain
#   yTrain=dataTrain$response
#   yTest=dataTest$response
#   
#   ## BRT
#   set.seed(1)
#   gbmOut=gbm(response ~ . ,data=dataTrain,
#              n.trees=hgrid$n.trees[row],
#              distribution="bernoulli",
#              shrinkage=hgrid$shrinkage[row],
#              interaction.depth=hgrid$interaction.depth[row],
#              n.minobsinnode=hgrid$n.minobsinnode[row],
#              cv.folds=5,class.stratify.cv=TRUE,
#              bag.fraction=0.5,train.fraction=1,
#              n.cores=1,
#              verbose=F)
#   
#   ## performance
#   par(mfrow=c(1,1),mar=c(4,4,1,1))
#   best.iter=gbm.perf(gbmOut,method="cv")
#   
#   ## predict with test data
#   preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")
#   
#   ## known
#   result=dataTest$response
#   
#   ## sensitiviy and specificity
#   sen=InformationValue::sensitivity(result,preds)
#   spec=InformationValue::specificity(result,preds)
#   
#   ## AUC on train
#   auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
#   
#   ## AUC on test
#   auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
#   
#   ## print
#   print(paste("hpar row ",row," done; test AUC is ",auc_test,sep=""))
#   
#   ## save outputs
#   return(list(best=best.iter,
#               trainAUC=auc_train,
#               testAUC=auc_test,
#               spec=spec,
#               sen=sen,
#               wrow=row))
# }
# 
# ## run the function
# hpars=lapply(1:nrow(hgrid),function(x) hfit(x))
# 
# ## export
# setwd("/Users/danielbecker/Desktop/prpred")
# saveRDS(hpars,"prpred_grid search 31923.rds")
# 
# ## get results
# hresults=data.frame(sapply(hpars,function(x) x$trainAUC),
#                     sapply(hpars,function(x) x$testAUC),
#                     sapply(hpars,function(x) x$spec),
#                     sapply(hpars,function(x) x$sen),
#                     sapply(hpars,function(x) x$wrow),
#                     sapply(hpars,function(x) x$best))
# names(hresults)=c("trainAUC","testAUC",
#                   "spec","sen","row","best")
# 
# ## combine and save
# search=merge(hresults,hgrid,by="row")
# 
# ## ggplot
# ggplot(search,aes(factor(interaction.depth),testAUC,group=factor(shrinkage),colour=factor(shrinkage)))+
#   geom_boxplot()+theme_bw()+
#   facet_grid(prop~n.minobsinnode)
# 
# ## aggregate
# smean=aggregate(cbind(trainAUC,testAUC,spec,sen)~n.trees+interaction.depth+shrinkage+n.minobsinnode+prop,
#                 data=search,mean)

## best pars = 90%, 4 ID, 0.001, nmin = 4

## function to use different data partitions
brt_part=function(seed,response){
  
  ## make new data
  ndata=set
  
  ## clean
  ndata$tip=NULL
  
  ## correct response
  ndata$response=ndata[response][,1]
  
  ## remove raw
  ndata$Pr=NULL
  
  ## use rsample to split 90/10
  set.seed(seed)
  split=initial_split(ndata,prop=0.90,strata="response")
  
  ## test and train
  dataTrain=training(split)
  dataTest=testing(split)
  
  ## yTest and yTrain
  yTrain=dataTrain$response
  yTest=dataTest$response

  ## BRT
  set.seed(1)
  gbmOut=gbm(response ~ . ,data=dataTrain,
             n.trees=6000,
             distribution="bernoulli",
             shrinkage=0.001,
             interaction.depth=4,
             n.minobsinnode=4,
             cv.folds=5,class.stratify.cv=TRUE,
             bag.fraction=0.5,train.fraction=1,
             n.cores=1,
             verbose=F)
  
  ## performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  best.iter=gbm.perf(gbmOut,method="cv")
  
  ## predict with test data
  preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")
  
  ## AUC on train
  auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
  
  ## AUC on test
  auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
  
  ## inner loop if yTest is all 0
  if(var(yTest)==0){
    
    perf=NA
  }else{
      
    ## ROC
    pr=prediction(preds,dataTest$response)
    perf=performance(pr,measure="tpr",x.measure="fpr")
    perf=data.frame(perf@x.values,perf@y.values)
    names(perf)=c("fpr","tpr")
      
    ## add seed
    perf$seed=seed
      
    }
  
  ## relative importance
  bars=summary(gbmOut,n.trees=best.iter,plotit=F)
  bars$rel.inf=round(bars$rel.inf,2)
  
  ## temporary dataset to fix citations at the mean
  temp=set
  temp$lcites=mean(temp$lcites,na.rm=T)
  
  ## predict with cites included
  preds=predict(gbmOut,set,n.trees=best.iter,type="response")
  pred_data=raw[c("tip","zero","Pr","English")]
  pred_data$pred=preds
  pred_data$type=response
  
  ## fix citations
  pred_data$pred_cites=predict(gbmOut,temp,n.trees=best.iter,type="response")
  
  ## sort
  pred_data=pred_data[order(pred_data$pred_cites,decreasing=T),]
  
  ## print
  print(paste("BRT ",seed," done; test AUC = ",round(auc_test,3),sep=""))
  
  ## save outputs
  return(list(mod=gbmOut,
              best=best.iter,
              trainAUC=auc_train,
              testAUC=auc_test,
              roc=perf,
              rinf=bars,
              predict=pred_data,
              traindata=dataTrain,
              testdata=dataTest,
              seed=seed))
}

## run function
smax=100
brts=lapply(1:smax,function(x) brt_part(seed=x,response="Pr"))

## mean test AUC
round(mean(sapply(brts,function(x) x$testAUC))*100,0)
round(se(sapply(brts,function(x) x$testAUC))*100,2)
## 80% accuracy +/- 1.74%

## relative importance
vinf=lapply(brts,function(x) x$rinf)
vinf=do.call(rbind,vinf)

## aggregate mean
vdata=data.frame(aggregate(rel.inf~var,data=vinf,mean),
                 aggregate(rel.inf~var,data=vinf,se)["rel.inf"])
names(vdata)=c("var","rel.inf","rse")
vdata=vdata[order(vdata$rel.inf,decreasing=T),]
rm(vinf)

## make rmin and rmax
vdata$rmin=vdata$rel.inf-vdata$rse
vdata$rmax=vdata$rel.inf+vdata$rse

## trim
vset=vdata
vset=vset[vset$rel.inf>1,]
vset=vset[order(vset$rel.inf),]
vset$var=factor(vset$var,levels=unique(vset$var))

## plot
ggplot(vset,aes(reorder(var,rel.inf,max),rel.inf))+
  geom_segment(aes(y=0,yend=rel.inf,
                   x=var,xend=var))+
  geom_point(size=1.5)+
  #scale_y_sqrt()+
  coord_flip()+
  theme_bw()+
  labs(x=NULL,
       y="mean feature importance")+
  theme(axis.title=element_text(size=10),
        axis.text=element_text(size=8))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))

## average predictions
apreds=lapply(brts,function(x) x$predict)
apreds=do.call(rbind,apreds)

## aggregate
apreds=data.frame(aggregate(pred~tip,data=apreds,mean))

## get basic info
dat=raw[c("tip","English","Pr","zero","Migration","urban","Family3","Order3","IUCN")]

## merge
apreds=merge(apreds,dat,by="tip")

## threshold
to=optimal.thresholds(data.frame(apreds[c('tip','Pr','pred')]),
                      threshold = 10001,
                      opt.methods = 10,
                      req.sens = 0.95,
                      na.rm = TRUE)

## binary
apreds$bin=ifelse(apreds$pred>to$pred,1,0)
table(apreds$bin)

## make known/unknown

## ridgeline
library(ggridges)
ggplot(apreds[!is.na(apreds$IUCN),],
       aes(pred,IUCN))+
  geom_density_ridges()+
  facet_wrap(~Pr)+
  scale_x_continuous(trans="log10")

