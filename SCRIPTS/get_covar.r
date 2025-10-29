# args <- c("DATA/geno/2024-04-19_IMAGEcomplete_CR_SCM_GGA6_Ref0_Alt1.RData",
#     "INTERMEDIATE/covar/covar.txt", "INTERMEDIATE/covar/covar_bfmap.csv")

args <- commandArgs(TRUE)
args

library(data.table)

load(args[1])
PED <- setDT(PED)
covar <- PED[, .(Patient_ID, Generation)]


write.table(cbind("IM", covar), args[2], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.csv(covar, args[3], quote = FALSE, row.names = FALSE)