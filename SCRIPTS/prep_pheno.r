# for gcta analysis the phenotype file: famID, IID, phenotypes columns

# args <-  c("DATA/geno/2024-04-19_IMAGEcomplete_CR_SCM_GGA6_Ref0_Alt1.RData", "DATA/pheno/Phenotypes_Final.csv", "DATA/pheno/BW32.txt", "BW32")
args <- commandArgs(TRUE)
args

library(data.table)
load(args[1])
pheno <- fread(args[2])

colnames(pheno)
colnames(PED)
colnames(PED)[5] = "kml"
trait = args[4]
meta = merge(PED[,c("Patient_ID","kml")], pheno[,c("kml",..trait)], by = "kml")
colnames(meta)
setDT(meta)
meta = meta[!which(is.na(meta[,3])),]
meta[, kml := Patient_ID]
write.table(meta, args[3], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
