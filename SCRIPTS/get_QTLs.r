# args <- c("RESULTS/gcta/mlma_impute/all/BW32_all.mlma", "BW32", "RESULTS/metaGWAS/ENs_all.txt",
# "RESULTS/BFMAP/BW32/QTLs_BW32_all.txt","11000000")

args <- commandArgs(TRUE)
args

library(data.table)
options(scipen = 999)

mlmap <- args[1]
pheno <- args[2]
mgwap <- args[3]
outxt <- args[4]
maxlength <- as.numeric(args[5])
testphenos <- c("BW32", "EW30", "EW70", "EN1", "EN13")

# QTL from multi-trait GWAS (EN1 - EN13)
ens <- paste0("EN",1:13)

print(pheno)

if (pheno %in% ens) {
    mlma <- fread(mgwap)
} else {
    mlma <- fread(mlmap)
}


bonf <- 0.05 / nrow(mlma)
sig <- mlma[p < bonf,]
chrs <- unique(sig$Chr)

st <- sig[, .SD[1], by = .(Chr)][,1:3]
ed <- sig[, .SD[.N], by = .(Chr)][,1:3]
tab <- st[ed, on = .(Chr = Chr)]
setnames(tab, "i.SNP", "edSNP")
setnames(tab, "i.bp", "edbp")
tab[, exdst := bp - 1e6]
tab[, exded := edbp + 1e6]

if (pheno %in% testphenos) {

    tab <- data.table(Chr = tab$Chr,
                exdst = c(tab$exdst, tab$exdst - 1e6, tab$exdst - 2e6),
                exded = c(tab$exded, tab$exded + 1e6, tab$exded + 2e6)
    )
}

tab[, diff := exded - exdst]


# for (c in chrs){
#     if (tab[Chr == c, diff] > maxlength){
#         print(paste0("chr ",c," interval > 10e6"))
#         sigc <- sig[Chr == c]
#         sigc[, window := floor(bp / 10e6)]
#         if (pheno == "EN2"){
#             sigc[window == 8, window := 9]
#         }
#         st <- sigc[, .SD[1], by = .(Chr, window)][,1:4]
#         ed <- sigc[, .SD[.N], by = .(Chr, window)][,1:4]
#         sted <- st[ed, on = .(Chr = Chr, window = window)]
#         setnames(sted, "i.SNP", "edSNP")
#         setnames(sted, "i.bp", "edbp")
#         sted[, exdst := bp - 1e6]
#         sted[, exded := edbp + 1e6]
#         sted[, diff := exded - exdst]
#         sted <- sted[, -c("window")]
        
#         tab <- tab[!Chr == c,]
#         tab <- rbind(tab, sted)
#     }
# }
tab <- tab[order(Chr, exdst)]

write.table(tab, outxt, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
