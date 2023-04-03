## Plasmodium relictum prediction (prpred)
## 02_trait merge
## danbeck@ou.edu
## last updated 4/2/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## libraries
library(plyr)
library(readxl)
library(ape)
library(phyloregion)

## load flat files
setwd("/Users/danielbecker/Desktop/prpred/MalAvi flat files")
data=read.csv("binary Plasmodium relictum status MalAvi.csv")
data$X=NULL

## fix Jetz
data$Jetz.species=ifelse(data$Jetz.species=="",NA,data$Jetz.species)

## if Jetz match, use that name
data$tip=ifelse(is.na(data$Jetz.species),data$species,data$Jetz.species)

## aggregate
set=aggregate(Pr~tip,data=data,sum)
set$Pr=ifelse(set$Pr>0,1,0)

## bind back with data
data=data[!duplicated(data$tip),]
data$Pr=NULL
data=merge(data,set,by="tip")
rm(set)

## fix tip
data$tip=gsub("_"," ",data$tip)

## load BirdTree taxonomy for backbone
setwd("~/OneDrive - University of Oklahoma/Becker Lab/Datasets/BirdTree")
pnames=read.csv("BLIOCPhyloMasterTax.csv")

## check mismatch
mis=setdiff(data$tip,pnames$Scientific) ## 80 mismatch

## manually edit remaining 80 mismatch
data$tip=revalue(data$tip,
                 c("Acritillas indica"="Iole indica",
                   "Actitis macularia"="Actitis macularius",
                   "Anisognathus flavinuchus"="Anisognathus somptuosus",
                   "Anser canagica"="Chen canagica",
                   "Anthropoides paradiseus"="Grus paradisea",
                   "Anthus cinnamomeus"="Anthus novaeseelandiae",
                   "Bubo scandiacus"="Bubo scandiaca",
                   "Cantorchilus longirostris"="Thryothorus longirostris",
                   "Cantorchilus superciliaris"="Thryothorus superciliaris",
                   "Cardellina canadensis"="Wilsonia canadensis",
                   "Chelidorhynx hypoxanthus"="Chelidorhynx hypoxantha",
                   "Chloris sinica"="Carduelis sinica",
                   "Chlorospingus flavopectus"="Chlorospingus ophthalmicus",
                   "Chroicocephalus ridibundus"="Larus ridibundus",
                   "Cichlherminia lherminieri"="Turdus lherminieri",
                   "Clanga pomarina"="Aquila pomarina",
                   "Corythornis madagascariensis"="Ceyx madagascariensis",
                   "Criniger pallidus"="Alophoixus pallidus",
                   "Cyanoloxia brissonii"="Cyanocompsa brissonii",
                   "Cyornis pallidipes"="Cyornis pallipes",
                   "Cypseloides rutilus"="Streptoprocne rutila",
                   "Daptrius americanus"="Ibycter americanus",
                   "Diglossopis glauca"="Diglossa glauca",
                   "Doryfera ludoviciae"="Doryfera ludovicae",
                   "Elaenia chilensis"="Elaenia albiceps",
                   "Erythrogenys erythrogenys"="Pomatorhinus erythrogenys",
                   "Euaegotheles insignis"="Aegotheles insignis",
                   "Falcipennis canadensis"="Dendragapus canadensis",
                   "Geokichla citrina"="Zoothera citrina",
                   "Geospizopsis unicolor"="Phrygilus unicolor",
                   "Griseotyrannus aurantioatrocristatus"="Empidonomus aurantioatrocristatus",
                   "Gymnoris pyrgita"="Petronia pyrgita",
                   "Hartlaubius auratus"="Saroglossa aurata",
                   "Henicorhina anachoreta"="Henicorhina leucophrys",
                   "Hesperiphona vespertina"="Coccothraustes vespertinus",
                   "Ispidina picta"="Ceyx pictus",
                   "Ketupa ketupa"="Ketupa ketupu",
                   "Kittacincla malabarica"="Copsychus malabaricus",
                   "Machlolophus spilonotus"="Parus spilonotus",
                   "Malurus assimilis"="Malurus lamberti",
                   "Melanitta americana"="Melanitta nigra",
                   "Melanitta stejnegeri"="Melanitta fusca",
                   "Microscelis amaurotis"="Ixos amaurotis",
                   "Mixornis gularis"="Macronous gularis",
                   "Momotus lessonii"="Momotus aequatorialis",
                   "Montecincla jerdoni"="Strophocincla cachinnans",
                   "Montecincla meridionalis"="Strophocincla fairbanki",
                   "Myiothlypis flaveola"="Basileuterus flaveolus",
                   "Nyctipolus nigrescens"="Caprimulgus nigrescens",
                   "Oedistoma iliolophus"="Toxorhamphus iliolophus",
                   "Origma murina"="Crateroscelis murina",
                   "Origma robusta"="Crateroscelis robusta",
                   "Parus minor"="Parus major",
                   "Peneothello sigillata"="Peneothello sigillatus",
                   "Pheugopedius eisenmanni"="Thryothorus eisenmanni",
                   "Philohydor lictor"="Pitangus lictor",
                   "Piculus rivolii"="Colaptes rivolii",
                   "Piculus rubiginosus"="Colaptes rubiginosus",
                   "Poecile sclateri"="Parus sclateri",
                   "Psittiparus gularis"="Paradoxornis gularis",
                   "Rhynchospiza stolzmanni"="Aimophila stolzmanni",
                   "Rupornis magnirostris"="Buteo magnirostris",
                   "Sciaphylax hemimelaena"="Myrmeciza hemimelaena",
                   "Seicercus xanthodryas"="Phylloscopus borealis",
                   "Setophaga citrina"="Wilsonia citrina",
                   "Sholicola ashambuensis"="Myiomela albiventris",
                   "Sholicola major"="Myiomela major",
                   "Silvicultrix frontalis"="Ochthoeca frontalis",
                   "Silvicultrix jelskii"="Ochthoeca jelskii",
                   "Spizaetus cirrhatus"="Nisaetus cirrhatus",
                   "Sporophila angolensis"="Oryzoborus angolensis",
                   "Stomiopera unicolor"="Lichenostomus unicolor",
                   "Syndactyla ucayalae"="Simoxenops ucayalae",
                   "Thamnophilus bernardi"="Sakesphorus bernardi",
                   "Tringa brevipes"="Heteroscelus brevipes",
                   "Trochalopteron cachinnans"="Strophocincla cachinnans",
                   "Trochalopteron fairbanki"="Strophocincla fairbanki",
                   "Turdus eunomus"="Turdus naumanni",
                   "Tyto furcata"="Tyto alba",
                   "Zoothera aurea"="Zoothera dauma"))

## recheck
mis=setdiff(data$tip,pnames$Scientific)
rm(mis)

## aggregate by tip
set=aggregate(Pr~tip,data=data,sum)
set$Pr=ifelse(set$Pr>0,1,0)

## bind back with data
data=data[!duplicated(data$tip),]
data$Pr=NULL
data=merge(data,set,by="tip")
rm(set)

## merge data with all birds
pnames$tip=pnames$Scientific
data=merge(pnames,data,by="tip",all=T)

## clean
rm(pnames)

## assign pseudoabsences
data$zero=ifelse(is.na(data$Pr),"pseudoabsence","true zero")
data$Pr=ifelse(is.na(data$Pr),0,data$Pr)

## load BirdTree consensus phylogeny
setwd("/Users/danielbecker/Desktop/prpred/BirdTree phylo")
tree=readRDS("BirdNet 2K tree consensus.rds")

## fix labels
tree$tip.label=gsub("_"," ",tree$tip.label)

## check
mis=setdiff(tree$tip.label,data$tip)
rm(mis)

## evolutionary distinctiveness
a=phyloregion::evol_distinct(tree,type="fair.proportion")
b=phyloregion::evol_distinct(tree,type="equal.splits")

## check correlation
cor(a,b)

## save
a=data.frame(ed_fair=a)
b=data.frame(ed_equal=b)

## names
a$tip=rownames(a)
b$tip=rownames(b)

## merge
edata=merge(a,b,by="tip")
rm(a,b)

## combine into data
data=merge(data,edata,by="tip")
rm(edata)

## load AVONET
setwd("/Users/danielbecker/OneDrive - University of Oklahoma/Becker Lab/Datasets/AVONET")
avonet=read_excel("AVONET Supplementary dataset 1.xlsx",4)

## check name
mis=setdiff(data$tip,avonet$Species3)

## fix label
avonet$tip=avonet$Species3

## merge
data=merge(data,avonet,by="tip",all.x=T)

## Merge in Gonzalez-Lagos et al. 2022 (urban traits)
setwd("/Users/danielbecker/OneDrive - University of Oklahoma/Becker Lab/Datasets/Gonzalez-Lagos et al 2021 Ecography")
gdata=read.csv("DataS1.csv")

## fix
gdata$tip=gsub("_"," ",gdata$animal)

## check name
mis=setdiff(gdata$tip,data$tip)

## merge
gdata$animal=NULL
gdata=gdata[c("tip","urban","humanDisturbed","invasion.potential","UTI","tolerance")]
data=merge(data,gdata,by="tip",all.x=T)
rm(gdata,mis,avonet)

## export Pr with traits
setwd("/Users/danielbecker/Desktop/prpred")
write.csv(data,"MalAvi with host traits.csv")