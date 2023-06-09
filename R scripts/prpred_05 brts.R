## Plasmodium relictum prediction (prpred)
## 05_brts
## danbeck@ou.edu
## last updated 4/6/2023

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

## tabulate when trimming by order vs family
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
set$status=NULL

## tally
ncol(set)-2 ## tip and Pr

## factors
set$Migration=factor(set$Migration)
set$Trophic.Level=factor(set$Trophic.Level)
set$Primary.Lifestyle=factor(set$Primary.Lifestyle)
set$urban=factor(set$urban)
set$humanDisturbed=factor(set$humanDisturbed)
set$strategy_3=NULL
set$distance_4=factor(set$distance_4)

## extract continuous vars
vars=names(which(unlist(lapply(set,class))!="factor"))
vars=vars[!vars%in%c("tip","Pr","status")]

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

## hyperparameter grid
hgrid=expand.grid(n.trees=5000,
                  interaction.depth=c(2,3,4),
                  shrinkage=c(0.1,0.01),
                  n.minobsinnode=c(4,5),
                  seed=seq(1,5,by=1),
                  prop=c(0.8,0.9))

## trees, depth, shrink, min, prop
hgrid$id=with(hgrid,paste(n.trees,interaction.depth,shrinkage,n.minobsinnode))

## sort by id then seed
hgrid=hgrid[order(hgrid$id,hgrid$seed),]

## now add rows
hgrid$row=1:nrow(hgrid)

## factor id
hgrid$id2=factor(as.numeric(factor(hgrid$id)))

## function to assess each hyperpar combination
hfit=function(row){

  ## make new data
  ndata=set
  ndata$tip=NULL

  ## correct response
  ndata$response=ndata$Pr
  ndata$Pr=NULL

  ## use rsample to split
  set.seed(hgrid$seed[row])
  split=initial_split(ndata,prop=hgrid$prop[row],strata="response")

  ## test and train
  dataTrain=training(split)
  dataTest=testing(split)

  ## yTest and yTrain
  yTrain=dataTrain$response
  yTest=dataTest$response

  ## BRT
  set.seed(1)
  gbmOut=gbm(response ~ . ,data=dataTrain,
             n.trees=hgrid$n.trees[row],
             distribution="bernoulli",
             shrinkage=hgrid$shrinkage[row],
             interaction.depth=hgrid$interaction.depth[row],
             n.minobsinnode=hgrid$n.minobsinnode[row],
             cv.folds=5,class.stratify.cv=TRUE,
             bag.fraction=0.5,train.fraction=1,
             n.cores=1,
             verbose=F)

  ## performance
  # par(mfrow=c(1,1),mar=c(4,4,1,1))
  # best.iter=gbm.perf(gbmOut,method="cv")
  best.iter=gbm.perf(gbmOut,method="cv",plot.it=F)

  ## predict with test data
  preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")

  ## known
  result=dataTest$response

  ## sensitiviy and specificity
  sen=InformationValue::sensitivity(result,preds)
  spec=InformationValue::specificity(result,preds)

  ## AUC on train
  auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))

  ## AUC on test
  auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))

  ## print
  print(paste("hpar row ",row," done; test AUC is ",auc_test,sep=""))

  ## save outputs
  return(list(best=best.iter,
              trainAUC=auc_train,
              testAUC=auc_test,
              spec=spec,
              sen=sen,
              wrow=row))
}

## run the function
hpars=lapply(1:nrow(hgrid),function(x) hfit(x))

## export
setwd("/Users/danielbecker/Desktop/prpred")
saveRDS(hpars,"prpred_grid search 42123.rds")

## get results
hresults=data.frame(sapply(hpars,function(x) x$trainAUC),
                    sapply(hpars,function(x) x$testAUC),
                    sapply(hpars,function(x) x$spec),
                    sapply(hpars,function(x) x$sen),
                    sapply(hpars,function(x) x$wrow),
                    sapply(hpars,function(x) x$best))
names(hresults)=c("trainAUC","testAUC",
                  "spec","sen","row","best")

## combine and save
search=merge(hresults,hgrid,by="row")

## export
setwd("/Users/danielbecker/Desktop/prpred")
write.csv(search,"prpred_grid search 42123.csv")

## ggplot
ggplot(search,aes(factor(interaction.depth),testAUC,colour=factor(shrinkage)))+
  geom_boxplot()+theme_bw()+
  facet_grid(prop~n.minobsinnode)

## linear model test
mod=lm(testAUC~factor(interaction.depth)+factor(shrinkage)+factor(prop)+factor(n.minobsinnode),data=search)
anova(mod)
summary(mod)
## shrinkage = 0.01, prop = 0.8, ID = 4

## aggregate
smean=aggregate(cbind(trainAUC,testAUC,spec,sen)~n.trees+interaction.depth+shrinkage+n.minobsinnode+prop,
                data=search,mean)

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
  split=initial_split(ndata,prop=0.80,strata="response")
  
  ## test and train
  dataTrain=training(split)
  dataTest=testing(split)
  
  ## yTest and yTrain
  yTrain=dataTrain$response
  yTest=dataTest$response

  ## BRT
  set.seed(1)
  gbmOut=gbm(response ~ . ,data=dataTrain,
             n.trees=5000,
             distribution="bernoulli",
             shrinkage=0.01,
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
    perf=ROCR::performance(pr,measure="tpr",x.measure="fpr")
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
smax=50
brts=lapply(1:smax,function(x) brt_part(seed=x,response="Pr"))

## mean test AUC
round(mean(sapply(brts,function(x) x$testAUC))*100,0)
round(se(sapply(brts,function(x) x$testAUC))*100,2)
## 90% accuracy +/- 0.3%

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

## pdp
detach("package:purrr", unload=TRUE)
library(pdp)

## function for compiling across BRTs for a given predictor, all else equal
pdp_agg=function(mod,feature){
  
  ## just the plot function
  pdep=plot(mod$mod,feature,
            return.grid=T,
            n.trees=mod$best,
            plot=F,
            continuous.resolution=200,
            type="response")
  
  ## add seed
  pdep$seed=unique(mod$roc$seed)
  
  ## save predictor
  pdep$predictor=pdep[feature][,1]
  
  ## order
  pdep=pdep[order(pdep$predictor),]
  
  ## get rank
  pdep$rank=1:nrow(pdep)
  
  ## save yhat
  pdep$yhat=pdep$y
  
  ## return
  return(pdep)
  
}

## function to plot
pdp_plot=function(bmods,feature){
  
  ## pdp_agg
  agg=do.call(rbind,lapply(bmods,function(x) pdp_agg(x,feature)))
  
  ## get class of the feature
  cl=class(set[feature][,1])
  
  ## if else based on type
  if(cl%in%c("numeric","integer")){
    
    ## get element-wise means
    x=with(agg,tapply(predictor,rank,mean))
    y=with(agg,tapply(yhat,rank,mean))
    
    ## save as mean
    pmean=data.frame(predictor=x,yhat=y)
    
    ## get yrange
    yrange=range(agg$yhat,pmean$yhat,na.rm=T)
    
    ## get histogram
    hi=hist(raw[feature][,1],breaks=30,plot=F)
    hi=with(hi,data.frame(breaks[1:(length(breaks)-1)],counts))
    names(hi)=c("mids","counts")
    
    ## ggplot it
    ggplot(agg,aes(predictor,yhat,group=seed))+
      
      ## add histogram
      geom_segment(data=hi,inherit.aes=F,
                   aes(x=mids,xend=mids,
                       y=yrange[1],yend=plotrix::rescale(counts,yrange)),
                   size=1,colour="grey",alpha=0.25)+
      
      ## add lines
      geom_line(size=1,alpha=0.25,colour="grey")+
      
      ## add mean
      geom_line(data=pmean,size=2,inherit.aes=F,
                aes(predictor,yhat))+
      
      ## theme
      theme_bw()+
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=7))+
      theme(axis.title.x=element_text(margin=margin(t=5,r=0,b=0,l=0)))+
      theme(axis.title.y=element_text(margin=margin(t=0,r=5,b=0,l=0)))+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      labs(x=feature,y="marginal effect")+
      scale_y_continuous(labels=scales::number_format(accuracy=0.01))
    
    ## end numeric
  }else{ ## factor-based plot
    
    ## get element-wise means
    y=with(agg,tapply(yhat,predictor,mean))
    
    ## save as mean
    #pmean=data.frame(predictor=x,yhat=y)
    pmean=data.frame(y)
    names(pmean)="yhat"
    pmean$predictor=rownames(pmean)
    rownames(pmean)=NULL
    
    ## make temp data
    temp=set
    temp$predictor=temp[feature][,1]
    
    ## do nothing
    agg=agg
    pmean=pmean
    temp=temp
    
    ## get yrange
    yrange=range(agg$yhat,pmean$yhat,na.rm=T)
    
    ## fix temp to yrange
    temp$yhat=ifelse(temp$Pr==1,max(yrange),min(yrange))
    
    ## ggplot with rug
    set.seed(1)
    ggplot(agg,aes(predictor,yhat,group=seed))+
      
      ## add individual BRTs
      geom_jitter(size=1,alpha=0.25,colour="grey",width=0.1)+
      
      ## add mean
      geom_point(data=pmean,size=2,inherit.aes=F,shape=15,
                 aes(predictor,yhat))+
      
      ## add rug
      geom_rug(data=temp,inherit.aes=F,
               aes(predictor,yhat),
               sides="b",position="jitter",
               colour="grey",alpha=0.25,
               na.rm=T)+
      
      ## theme
      theme_bw()+
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=7))+
      theme(axis.title.x=element_text(margin=margin(t=5,r=0,b=0,l=0)))+
      theme(axis.title.y=element_text(margin=margin(t=0,r=5,b=0,l=0)))+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      labs(x=feature,y="marginal effect")+
      scale_y_continuous(limits=c(yrange[1]-0.01,yrange[2]+0.01),
                         labels=scales::number_format(accuracy=0.01))
    
  }
  
}

## visualize
pdp_plot(brts,vdata$var[1])
pdp_plot(brts,vdata$var[2])
pdp_plot(brts,vdata$var[3])
pdp_plot(brts,vdata$var[4])
pdp_plot(brts,vdata$var[5])
pdp_plot(brts,vdata$var[6])
pdp_plot(brts,vdata$var[7])
pdp_plot(brts,vdata$var[8])
pdp_plot(brts,vdata$var[9])
