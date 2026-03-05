# args = c("DATA/pheno/Phenotypes_Final.csv", "DATA/geno/2024-04-19_IMAGEcomplete_CR_SCM_GGA6_Ref0_Alt1.RData")
args = commandArgs(TRUE)
args
library(data.table)
pheno = fread(args[1])
load(args[2])
tt = c('EN1','EN2','EN3','EN4','EN5','EN6','EN7','EN8','EN10','EN11','EN12','EN13','BW32','EW30','EW40','EW50','EW70')


# --phenotype_file: bfmap requires csv, first row header, first col ID, other cols phenos, missing leave empty
PED = setDT(PED)
PED$Patient_KM = as.integer(PED$Patient_KM)
pp = pheno[PED[,c("Patient_ID","Patient_KM")], on=.(kml = Patient_KM), nomatch=NULL]
dim(pp)[1] == sum(pheno$kml %in% PED$Patient_KM)
cols = c("Patient_ID", tt)
pp = pp[, ..cols]
colnames(pp)

write.csv(pp, "tmp_pheno.csv", row.names=FALSE, quote=FALSE, na="")

