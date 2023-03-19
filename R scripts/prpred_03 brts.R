## Plasmodium relictum prediction (prpred)
## 03_brts
## danbeck@ou.edu
## last updated 3/18/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## libraries
library(gbm)
library(fastDummies)
library(rsample)
library(ROCR)

## load in MalAvi with traits
setwd("/Users/danielbecker/Desktop/prpred")
data=read.csv("MalAvi with host traits.csv")

## get orders of positive birds
pos=data[data$Pr==1,]
porder=unique(pos$Order3) ## AVONET

## trim
data=data[data$Order3%in%porder,]

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
data$Species.Status
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

## save raw
raw=data

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

## mode function
mode.prop <- function(x) {
  ux <- unique(x[is.na(x)==FALSE])
  tab <- tabulate(match(na.omit(x), ux))
  max(tab)/length(x[is.na(x)==FALSE])
}

## assess variation across columns
vars=data.frame(apply(data,2,function(x) mode.prop(x)),
                apply(data,2,function(x) length(unique(x))))

## get names
vars$variables=rownames(vars)
names(vars)=c("var","uniq","column")
vars$var=round(vars$var,2)

## if homogenous (100%)
vars$keep=ifelse(vars$var==1,"cut","keep")
vars=vars[order(vars$keep),]

## trim
keeps=vars[-which(vars$keep=="cut"),]$column

## drop if no variation
data=data[keeps]
rm(keeps,vars)

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
rm(mval,keeps,porder,i,mode.prop)

## set
set=data
set$TipLabel=NULL
set$Hackett_FineClades=NULL
set$Hackett_CoarseClades=NULL
set$English=NULL
set$BLFamilyLatin=NULL
set$Species3=NULL
set$zero=NULL

## factors
set$Habitat=factor(set$Habitat)
set$Migration=factor(set$Migration)
set$Trophic.Level=factor(set$Trophic.Level)
set$Trophic.Niche=factor(set$Trophic.Niche)
set$Primary.Lifestyle=factor(set$Primary.Lifestyle)
set$urban=factor(set$urban)
set$humanDisturbed=factor(set$humanDisturbed)
set$IUCN=factor(set$IUCN)

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
  
  ## use rsample to split 80/20
  set.seed(seed)
  split=initial_split(ndata,prop=0.8,strata="response")
  
  ## test and train
  dataTrain=training(split)
  dataTest=testing(split)
  
  ## yTest and yTrain
  yTrain=dataTrain$response
  yTest=dataTest$response

  ## BRT
  set.seed(1)
  gbmOut=gbm(response ~ . ,data=dataTrain,
             n.trees=1000,
             distribution="bernoulli",
             shrinkage=0.001,
             interaction.depth=3,
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
  
  ## predict
  preds=predict(gbmOut,set,n.trees=best.iter,type="response")
  pred_data=raw[c("tip","zero","Pr","English")]
  pred_data$pred=preds
  pred_data$type=response
  
  ## sort
  pred_data=pred_data[order(pred_data$pred,decreasing=T),]
  
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
smax=10
brts=lapply(1:smax,function(x) brt_part(seed=x,response="Pr"))
