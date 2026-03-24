
# args <- c("RESULTS/supdata/supdata.xlsx", 
# "RESULTS/BFMAPrecalc/BWEWs/credsets_BWEWs_all.txt",
# "RESULTS/BFMAPrecalc/ENs/credsets_ENs_all.txt",
# "RESULTS/FINEMAP/EN1/all/sss/credsets.txt",
# "RESULTS/FINEMAP/EN13/all/sss/credsets.txt",
# "RESULTS/FINEMAP/BW32/all/sss/credsets.txt")


args <- commandArgs(TRUE)
args

library(data.table)
library(writexl)

outp <- args[1]
filesp <- args[2:length(args)]
# porder <- c('EN1','EN2','EN3','EN4','EN5','EN6','EN7','EN8','EN10','EN11','EN12','EN13','BW32','EW30','EW40','EW50','EW70')

#
filesp_bfmap_en <- filesp[grep("BFMAPrecalc/ENs", filesp)]
filesp_bfmap_w <- filesp[grep("BFMAPrecalc/BWEWs", filesp)]

filesp_finemap <- filesp[grep("FINEMAP", filesp)]
filesp_finemap_en <- filesp_finemap[grep("EN", filesp_finemap)]
filesp_finemap_w <- setdiff(filesp_finemap, filesp_finemap_en)

# combine bfmap
res11 <- rbindlist(lapply(filesp_bfmap_en, function(p) {

    dat <- fread(p)
}))

res12 <- rbindlist(lapply(filesp_bfmap_w, function(p) {

    dat <- fread(p)
}))

# combine finemap
res21 <- rbindlist(lapply(filesp_finemap_en, function(p) {

    dat <- fread(p)
    dat <- dat[log10bf >= 2 & prob > 0][, pheno := strsplit(p, "/")[[1]][3]]
}))

res22 <- rbindlist(lapply(filesp_finemap_w, function(p) {

    dat <- fread(p)
    dat <- dat[log10bf >= 2 & prob > 0][, pheno := strsplit(p, "/")[[1]][3]]
}))


#
out <- setNames(list(res11, res12, res21, res22), 
                c("BFMAP_ENs", "BFMAP_BWEWs", "FINEMAP_ENs", "FINEMAP_BWEWs"))
write_xlsx(out, path = outp)
