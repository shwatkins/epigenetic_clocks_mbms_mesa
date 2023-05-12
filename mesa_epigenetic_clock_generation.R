# epigenetic clocks creation
# MESA
# written by Sarah Watkins
# September 2021

# here we calculate 9 of the epigenetic clocks
# then we combine them all together in a df
# then we add in batch vars (slide + row)
# then we add in cell counts
# this is the dataset to send to Harvard for the epigenetic clocks
# analyses.

# Set directories:

project.dir <- "/project_dir"

data.dir <- file.path(project.dir, "data_dir")
output.dir <- file.path(project.dir,"output_dir")
clocks.dir <- file.path(project.dir, "clocks_dir")

# Creates directory (hash out once run)
#dir.create(output.dir)
output.dir

## http://github.com/perishky/meffil
library(meffil)
options(mc.cores=10)

## http://github.com/perishky/eval.save
library(eval.save)
eval.save.dir(output.dir)

library(readxl)
library(meffonym)

# load in data
# we need meth and we need age
# see dataset generation for creation of mesa_selected_variables.csv (phenotype data)
pheno <- read.csv(paste0(data.dir,"/phenotype/mesa_selected_variables.csv"))
# see conversion from M values to betas script for meth_betas.rds
meth <- readRDS(paste0(data.dir,"/omic/meth_betas.rds"))
meth <- as.data.frame(meth)
# add new clocks

#####################

# epiToc (385 CpG sites). "The average DNAm over these 385 CpG sites"
#meffonym.add.model("age.something", variables, coefficients, description)
epiToc <- read_xls(paste0(clocks.dir,"/13059_2016_1064_MOESM2_ESM.xls"))
epiToc <- epiToc[[1]][-1]
# not all of the epitoc cpgs may be present in our daatset:
present_epitoc_cpgs <- epiToc[epiToc%in% rownames(meth)]
length(present_epitoc_cpgs)
# 385
meth_epiToc <- meth[present_epitoc_cpgs,]
dim(meth_epiToc)
sum.cpgs <- colSums(meth_epiToc)
length(sum.cpgs)
epiToc.clock <- sum.cpgs/length(epiToc)
length(epiToc.clock)
head(epiToc.clock)
epiToc.clock <- as.data.frame(epiToc.clock)

######################

# Zhang mortality (10 CpG sites)
# We can have a continuous or 0-10 risk score
# "We .. used the fourth quartile value of cg08362785 and first 
# quartile values of other nine CpGs as the cutoff points, to define 
# aberrant methylation for each CpG (the exact cutoff points are 
# listed in Supplementary Table 4)"
zhang_cpgs <- c("cg01612140", "cg05575921", "cg06126421", "cg08362785",
                "cg10321156", "cg14975410", "cg19572487", "cg23665802",
                "cg24704287", "cg25983901")
zhang_present_on_epic <- zhang_cpgs[zhang_cpgs %in% rownames(meth)]
length(zhang_present_on_epic)
# 10
####### or make the continuous score
zhang_coefficients <- c(-0.38253, -0.92224, -1.70129, 2.71749,
                        -0.02073, -0.04156, -0.28069, -0.89440,
                        -2.98637, -1.80325)

meffonym.add.model("zhang.mortality", c("intercept", zhang_cpgs), c(0, zhang_coefficients), c("Zhang et al 2017 mortality clock from PMID 28303888. Clock details in supplement of the paper."))

########################

# MiAge (268 CpG sites)
miage <- read.csv(paste0(clocks.dir,"/2017EPI0166R-s02.csv"), stringsAsFactors=F)
# miage has functions and a script to get the mitotic age estimates
miage_in_epic <- miage$CpG_site_ID[miage$CpG_site_ID %in% rownames(meth)]
length(miage_in_epic)
#268 

source(paste0(clocks.dir,"/miage_function_library.r")) ### library of all functions used in the calculation of mitotic ages                                                                                 
clocksitesID=as.matrix(read.csv(paste0(clocks.dir,"/miage_Additional_File1.csv"),header=T))[,1]   ### epigenetic clock CpG sites
beta=  meth[ match(clocksitesID,rownames(meth)),]##select clock CpG sites 

load(paste0(clocks.dir,"/miage_site_specific_parameters.Rdata")) #### site-specific parameters

b=methyl.age[[1]];c=methyl.age[[2]];d=methyl.age[[3]]

miage.clock=mitotic.age(beta,b,c,d) ### estimated mitotic age
names(miage.clock)=colnames(beta)
head(miage.clock)
miage.clock <- as.data.frame(miage.clock)

#########################

# Levine/PhenoAge (513 sites)
levine <- read.csv(paste0(clocks.dir,"/aging-10-101414-s002.csv"),stringsAsFactors=F)
levine_in_epic <- levine$CpG[levine$CpG %in% rownames(meth)]
length(levine_in_epic)
# 513
levine_cpgs <- as.character(levine$CpG)
levine_weights <- levine$Weight
class(levine_weights)

meffonym.add.model("phenoage", variables=levine_cpgs, coefficients=levine_weights, c("PhenoAge (Levine) clock"))

meffonym.models()

#########################

# Dunedin
## devtools::install_github("danbelsky/DunedinPoAm38")
library("DunedinPoAm38")
dunedin.sites <- getRequiredProbes()$DunedinPoAm_38
length(dunedin.sites)
# 46
dunedin_in_epic <- dunedin.sites[dunedin.sites %in% rownames(meth)]
length(dunedin_in_epic)
# 46
# rows=cpgs
dunedin_age <- PoAmProjector(meth)
dunedin_age <- data.frame(dunedin_age[[1]])
names(dunedin_age) <- c("dunedin_age")
head(dunedin_age)

########################

# DNAmTL (140 CpG sites)
# load the sites and coefficients
dnamtl <- read.csv(paste0(clocks.dir,"/dnamtl_coefficients.csv"),stringsAsFactors=F)
colnames(dnamtl) <- c("cpg","weight")
head(dnamtl)
dnamtl_cpgs <- as.character(dnamtl$cpg)
length(dnamtl_cpgs)
dnamtl_in_epic <- dnamtl_cpgs[dnamtl_cpgs %in% rownames(meth)]
length(dnamtl_in_epic)
# 140
dnamtl_weights <- dnamtl$weight
class(dnamtl_weights)

meffonym.add.model("dnamtl", variables=dnamtl_cpgs, coefficients=dnamtl_weights, c("DNAmTL clock"))


######################

# Zhang (514 sites)
# script taken from github page

############# 2. data loading and QC ##################
print("1. Data loading and QC")

meth.t <- as.data.frame(t(meth))

addna<-function(methy){
  methy[is.na(methy)]<-mean(methy,na.rm=T)
  return(methy)
}

print("1.2 Replacing missing values with mean value")
#dataNona<-apply(data,2,function(x) addna(x))   ###############  replace the NA with mean value for each probe 

#dataNona<- dataNona[,apply(dataNona,2,function(x) sum(is.na(x)))!=nrow(dataNona)] ############ remove the probe when it has NA across all individuals
#print(paste0(ncol(data) - ncol(dataNona)," probe(s) is(are) removed since it has (they have) NA across all individuals"))

print("1.3 Standardizing")
meth.t<- apply(meth.t,1,scale)        ############### standardize the DNA methylation within each individual, remove the mean and divided by the SD of each individual     Probe * IND
rownames(meth.t)<-rownames(meth)

############# 3. get the coefficients of each probe from Elastic Net/BLUP method, !!!!WE HAVE TWO PREDICTORS!!!#############

print("2. Loading predictors")
read.table(paste0(clocks.dir,"/en.coef"),stringsAsFactor=F,header=T)->encoef

zhang_on_epic <- encoef$probe[encoef$probe %in% rownames(meth)]
length(zhang_on_epic)
#514
en_int<-encoef[1,2]

encoef<-encoef[-1,]

rownames(encoef)<-encoef$probe

############# 4. get common probes between predictors and data ##############
print("3. Checking misssing probes")

encomm<- intersect(rownames(encoef),rownames(meth.t))

endiff<- nrow(encoef) - length(encomm)

print(paste0(endiff," probe(s) in Elastic Net predictor is(are) not in the data"))

############# 5. extract the common probes and do age prediction ###############
print("4. Predicting")

encoef<-encoef[encomm,]

encoef$coef%*%meth.t[encomm,]+en_int->enpred

############# 6. Save the predicted result ###########

zhang_clock <- as.data.frame(t(enpred))
colnames(zhang_clock) <- c("zhang_clock")
head(zhang_clock)
summary(zhang_clock$zhang_clock)

# get Horvath and Hannum cpg numbers

horvath <- read.csv(paste0(clocks.dir,"/13059_2013_3156_MOESM3_ESM.csv"), stringsAsFactors=F, skip=2)
horvath <- horvath$CpGmarker[-1]
horvath_450k <- horvath[horvath %in% rownames(meth)]
length(horvath_450k)

hannum <- read_xlsx(paste0(clocks.dir,"/NIHMS418935-supplement-02.xlsx"))
hannum <- hannum$Marker
hannum_450k <- hannum[hannum %in% rownames(meth)]
length(hannum_450k)

## running meffonym clocks
rownames(pheno) <- pheno$idno
identical(rownames(pheno),colnames(meth))
meffonym.models()
meth.m <- as.matrix(meth)

meffonym_clocks <- cbind(age=pheno$age5c, ## check age colname
                         sapply(c("phenoage","zhang.mortality", "dnamtl","age.hannum","age.horvath"), function(model) {
                           meffonym.score(meth.m, model)$score
                         }))
head(meffonym_clocks)
meffonym_clocks <- as.data.frame(meffonym_clocks)
cor(meffonym_clocks)
summary(meffonym_clocks$phenoage)
summary(meffonym_clocks$zhang.mortality)
summary(meffonym_clocks$dnamtl)
rownames(meffonym_clocks) <- colnames(meth)
head(meffonym_clocks)

# add the other clocks
head(epiToc.clock)
head(miage.clock)
head(zhang_clock)
head(dunedin_age)

# epiToc.clock.450k, miage.clock.450k, dunedin_age.450k, zhang_clock.450k
all_clocks_mesa <- merge(meffonym_clocks, epiToc.clock, by="row.names")
head(all_clocks_mesa)
rownames(all_clocks_mesa) <- all_clocks_mesa$Row.names
head(all_clocks_mesa)
all_clocks_mesa <- all_clocks_mesa[,-1]
all_clocks_mesa <- merge(all_clocks_mesa, miage.clock, by="row.names")
rownames(all_clocks_mesa) <- all_clocks_mesa$Row.names
head(all_clocks_mesa)
all_clocks_mesa <- all_clocks_mesa[,-1]

all_clocks_mesa <- merge(all_clocks_mesa, dunedin_age, by="row.names")
rownames(all_clocks_mesa) <- all_clocks_mesa$Row.names
head(all_clocks_mesa)
all_clocks_mesa <- all_clocks_mesa[,-1]

all_clocks_mesa <- merge(all_clocks_mesa, zhang_clock, by="row.names")
rownames(all_clocks_mesa) <- all_clocks_mesa$Row.names
head(all_clocks_mesa)
all_clocks_mesa <- all_clocks_mesa[,-1]

head(all_clocks_mesa)
clock_seq <- colnames(all_clocks_mesa)
for(i in 1:ncol(all_clocks_mesa)){
  print(i)
  print(clock_seq[i])
  print(summary(all_clocks_mesa[,i]))
  print(sd(all_clocks_mesa[,i]))
}

save(all_clocks_mesa,file=paste0(output.dir,"/all_clocks_mesa.Rdata"))
cor(all_clocks_mesa)


# finally, load in batch info, add that to the clocks df,
# and add cell counts too

# load batch info
library(data.table)
batch <- fread(paste0(data.dir,"/omic/MESA_Epi_METH_idmap.txt"))
batch <- as.data.frame(batch)
rownames(batch) <- batch$idno
batch <- batch[,c("chip_meth", "pos_meth", "site")]
names(batch) <- c("slide", "row", "site")
head(batch)
batch$slide <- as.factor(as.character(batch$slide))
head(batch)
sum(is.na(batch))
identical(rownames(batch), colnames(meth)) 
identical(rownames(batch), rownames(all_clocks_mesa))
batch <- batch[,c(1:2)]
head(batch)
all_clocks_mesa <- merge(all_clocks_mesa,batch,by="row.names")
head(all_clocks_mesa)
rownames(all_clocks_mesa) <- all_clocks_mesa$Row.names
head(all_clocks_mesa)
all_clocks_mesa <- all_clocks_mesa[,-1]

# create cell counts
cellcounts<-meffil.estimate.cell.counts.from.betas(meth.m,cell.type.reference="blood gse35069 complete",verbose=T)
cellcounts<-data.frame(IID=row.names(cellcounts),cellcounts)
cellcounts <- cellcounts[,-1]
# make sure ID are row names!
head(cellcounts)
dim(cellcounts)
class(cellcounts$Bcell)
identical(rownames(all_clocks_mesa), rownames(cellcounts))
all_clocks_mesa <- merge(all_clocks_mesa,cellcounts,by="row.names")
head(all_clocks_mesa)
rownames(all_clocks_mesa) <- all_clocks_mesa$Row.names
head(all_clocks_mesa)
names(all_clocks_mesa)[names(all_clocks_mesa) == "Row.names"] <- "mesa_id"
head(all_clocks_mesa)

save(all_clocks_mesa,file=paste0(output.dir,"/all_clocks_mesa.Rdata"))
write.csv(all_clocks_mesa,file = paste0(output.dir,"/all_clocks_mesa.csv"),quote = F,row.names = T)
