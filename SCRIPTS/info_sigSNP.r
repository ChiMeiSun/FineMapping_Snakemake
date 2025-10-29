# args = c({input.mlma} {output.ss})
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
    
    print("min, max")
    sig_m[,min(bp), by=Chr][order(Chr)]
    sig_m[,max(bp), by=Chr][order(Chr)]
}

######
# folder <- c("mlma","fastGWA")
# gen <- c('all','offspring','FounderBC1','BC1','FounderBC2','BC2','IC')
# pheno <- paste0("EN",c(1:8,10:13))
# # pheno <- c("BW32","EW30","EW40","EW50","EW70")
# grid <- expand.grid(
#   folder = folder,
#   gen = gen,
#   pheno = pheno,
#   stringsAsFactors = FALSE
# )
# grid$file = sprintf("RESULTS/gcta/sigSNP/%s/%s/%s_%s.txt", grid$folder, grid$gen, grid$pheno, grid$gen)
# setDT(grid)
# grid[, sigChrs := vapply(file, function(f){
#   if (file.exists(f)) {
#     dt <- fread(f)
#     if ("Chr" %in% names(dt)) {
#       paste(sort(unique(dt$Chr)), collapse = ",")
#     } else {
#       NA_character_
#     }
#   } else {
#     NA_character_
#   }
# }, FUN.VALUE = character(1))]

# grid[, file := NULL]
# dat <- grid[!is.na(sigChrs)]
# dat <- dcast(dat, gen + pheno ~ folder, value.var = "sigChrs")
# dat[order(pheno)]
# dat[!is.na(fastGWA), .N, by = gen][order(-N)]
# dat[!is.na(mlma), .N, by = gen][order(-N)]