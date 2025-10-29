# args <- c("RESULTS/gcta/fastGWA/all/BW32_all.fastGWA", "RESULTS/metaGWAS/ENs_all.txt",
# "RESULTS/BFMAP/BW32/QTLs_BW32_all.txt","11000000")

args <- commandArgs(TRUE)
args

library(data.table)
options(scipen = 999)

# QTL from multi-trait GWAS (EN3 - EN13)
ens <- paste0("EN",1:13)
pheno <- unlist(strsplit(args[1],"[/_]"))[5]
print(pheno)

if (pheno %in% ens) {
    fgwa <- fread(args[2])
    colnames(fgwa) <- c("CHR","SNP","POS","A1","A2","FREQ","BETA","SE","P")
} else {
    fgwa <- fread(args[1])
}

bonf <- 0.05 / nrow(fgwa)
sig <- fgwa[P < bonf,]
chrs <- unique(sig$CHR)

st <- sig[, .SD[1], by = .(CHR)][,1:3]
ed <- sig[, .SD[.N], by = .(CHR)][,1:3]
tab <- st[ed, on = .(CHR = CHR)]
setnames(tab, "i.SNP", "edSNP")
setnames(tab, "i.POS", "edPOS")
tab[, exdst := POS - 1e6]
tab[, exded := edPOS + 1e6]
tab[, diff := exded - exdst]

maxlength <- as.numeric(args[4])
for (c in chrs){
    if (tab[CHR == c, diff] > maxlength){
        print(paste0("chr ",c," interval > 10e6"))
        sigc <- sig[CHR == c]
        sigc[, window := floor(POS / 10e6)]
        if (pheno == "EN2"){
            sigc[window == 8, window := 9]
        }
        st <- sigc[, .SD[1], by = .(CHR, window)][,1:4]
        ed <- sigc[, .SD[.N], by = .(CHR, window)][,1:4]
        sted <- st[ed, on = .(CHR = CHR, window = window)]
        setnames(sted, "i.SNP", "edSNP")
        setnames(sted, "i.POS", "edPOS")
        sted[, exdst := POS - 1e6]
        sted[, exded := edPOS + 1e6]
        sted[, diff := exded - exdst]
        sted <- sted[, -c("window")]
        
        tab <- tab[!CHR == c,]
        tab <- rbind(tab, sted)
    }
}
tab <- tab[order(CHR)]



write.table(tab, args[3], quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")