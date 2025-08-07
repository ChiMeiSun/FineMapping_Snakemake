# for gcta analysis the phenotype file: famID, IID, phenotypes columns

# args <-  {input.rdata} {input.pheno} {output.out} {wildcards.pheno}
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
meta[, kml := substr(Patient_ID, 1, 2)]
sum(is.na(meta[,3]))
meta = meta[!which(is.na(meta[,3])),]
dim(meta)

write.table(meta, args[3], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
