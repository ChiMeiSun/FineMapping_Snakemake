# args = c({input.mlma} {output.ss} {input.vcf})
args <- commandArgs(TRUE)
library(data.table)
m = fread(args[1], header=TRUE)
dim(m)

# methods: c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
# "BH" same as "FDR" -> Benjamini-Hochberg
# "BY" -> Benjamini-Yekutieli

# adjust p-value
# m$adjp = (p.adjust(m$p, method="bonferroni"))
# sig_m = m[m$adjp < 0.05,]
# sig_m = sig_m[order(sig_m$adjp),]

if (length(grep("fastGWA", args[1]) >1)){
    setnames(m, "P", "p")
    setnames(m, "CHR", "Chr")
    setnames(m, "POS", "bp")
}

bonf = 0.05/nrow(m)
sig_m = m[p < bonf,]
sig_m = sig_m[order(p),]
dim(sig_m)

if (nrow(sig_m) == 0) {
} else {
    write.table(sig_m, args[2], quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}

print("min, max")
sig_m[,min(bp), by=Chr][order(Chr)]
sig_m[,max(bp), by=Chr][order(Chr)]
