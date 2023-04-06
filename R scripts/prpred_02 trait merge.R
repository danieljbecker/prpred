## Plasmodium relictum prediction (prpred)
## 02_trait merge
## danbeck@ou.edu
## last updated 4/6/2023

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
edge=read.csv("MalAvi edgelist_cleaned.csv")
edge$X=NULL

## fix Jetz
edge$Jetz.species=ifelse(edge$Jetz.species=="",NA,edge$Jetz.species)

## if Jetz match, use that name
edge$tip=ifelse(is.na(edge$Jetz.species),edge$species,edge$Jetz.species)

## make new id
edge$id=paste(edge$tip,edge$Lineage_Name,sep="_")

## aggregate
set=aggregate(record~id,data=edge,sum)
edge=edge[!duplicated(edge$id),]
edge$record=NULL
edge=merge(set,edge,by="id",all.x=T)
rm(set)

## fix tip
edge$tip=gsub("_"," ",edge$tip)

## load BirdTree taxonomy for backbone
setwd("~/OneDrive - University of Oklahoma/Becker Lab/Datasets/BirdTree")
pnames=read.csv("BLIOCPhyloMasterTax.csv")

## check mismatch
mis=setdiff(edge$tip,pnames$Scientific) 

## manually edit mismatch
edge$tip=revalue(edge$tip,
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
mis=setdiff(edge$tip,pnames$Scientific)
rm(mis)

## make new id
edge$id=paste(edge$tip,edge$Lineage_Name,sep="_")

## aggregate
set=aggregate(record~id,data=edge,sum)
edge=edge[!duplicated(edge$id),]
edge$record=NULL
edge=merge(set,edge,by="id",all.x=T)
rm(set)

## export MalAvi edgelist harmonized
setwd("/Users/danielbecker/Desktop/prpred/MalAvi flat files")
write.csv(edge,"MalAvi edgelist_harmonized.csv")

## simplify to Pr records
pr_edge=edge[edge$Pr==1,]
write.csv(pr_edge,"MalAvi edgelist_harmonized_Pr.csv")

## simplify to known 0/1s
data=aggregate(Pr~tip,data=edge,sum)
data$Pr=ifelse(data$Pr>0,1,0)

## merge metadata
meta=edge[!duplicated(edge$tip),]
meta=meta[c("tip","order","family","genus","Jetz.species")]
data=merge(data,meta,by="tip")
rm(meta)

## merge data with all birds
pnames$tip=pnames$Scientific
data=merge(pnames,data,by="tip",all=T)

## clean
rm(pnames)

## assign pseudoabsences
data$zero=ifelse(is.na(data$Pr),"pseudoabsence","true zero")
data$Pr=ifelse(is.na(data$Pr),0,data$Pr)

## status
data$status=ifelse(data$zero=="true zero" & data$Pr==1,"positive",
                   ifelse(data$zero=="true zero" & data$Pr==0,"true zero",
                          "pseudoabsence"))

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

## combine with additional migration traits
setwd("/Users/danielbecker/OneDrive - University of Oklahoma/Becker Lab/Datasets/Dufour et al 2019 Journal of Biogeography")
dufor=read_xlsx("jbi13700-sup-0002-datas1.xlsx",1)

## tip
dufor$tip=gsub("_"," ",dufor$Sp.Scien.jetz)
mis=setdiff(data$tip,dufor$tip)
rm(mis)

## merge
data=merge(data,dufor[c("tip","strategy_3","distance_4","distance_quanti_RES0")],
           by="tip",all.x=T)
rm(dufor)

## export Pr with traits
setwd("/Users/danielbecker/Desktop/prpred")
write.csv(data,"MalAvi with host traits.csv")