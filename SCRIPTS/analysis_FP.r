# args <- c("RESULTS/analysis/BFMAPvsFINEMAP.txt",

# "RESULTS/BFMAP/BW32/FS/1:167395369-174490482_BW32_all.csv", 
# "RESULTS/BFMAP/BW32/FS/1:168395369-173490482_BW32_all.csv", 
# "RESULTS/BFMAP/BW32/FS/1:169395369-172490482_BW32_all.csv", #...   
# "RESULTS/BFMAP/BW32/FS/4:71674618-78718096_BW32_all.csv", 
# "RESULTS/BFMAP/BW32/FS/4:72674618-77718096_BW32_all.csv", 
# "RESULTS/BFMAP/BW32/FS/4:73674618-76718096_BW32_all.csv", 

# "RESULTS/BFMAP/EW70/FS/4:73699997-78296913_EW70_all.csv", 
# "RESULTS/BFMAP/EW70/FS/4:74699997-77296913_EW70_all.csv", 
# "RESULTS/BFMAP/EW70/FS/4:75699997-76296913_EW70_all.csv", 

# "RESULTS/FINEMAP/BW32/all/sss/1:167395369-174490482/data.snp", 
# "RESULTS/FINEMAP/BW32/all/sss/1:168395369-173490482/data.snp", 
# "RESULTS/FINEMAP/BW32/all/sss/1:169395369-172490482/data.snp",
# "RESULTS/FINEMAP/BW32/all/sss/4:71674618-78718096/data.snp",
# "RESULTS/FINEMAP/BW32/all/sss/4:72674618-77718096/data.snp",
# "RESULTS/FINEMAP/BW32/all/sss/4:73674618-76718096/data.snp",

# "RESULTS/FINEMAP/EW70/all/sss/4:73699997-78296913/data.snp",
# "RESULTS/FINEMAP/EW70/all/sss/4:74699997-77296913/data.snp",
# "RESULTS/FINEMAP/EW70/all/sss/4:75699997-76296913/data.snp"

# )

args <- commandArgs(TRUE)
args

library(data.table)

#
outxt <- args[1]
filesp <- args[2:length(args)]
length(filesp)

# functions
extract_info <- function(p) {

    if (length(grep("BFMAP", p)) > 0 ) {
        data.table(
            file = p,
            data = sub("^[^_]*_(.*)\\..*", "\\1", p), 
            region = sub(".*FS/([^_]+)_.*", "\\1", p),
            chr = sub(".*FS/(.*):.*", "\\1", p),
            method = "BFMAP-FS")

    } else if (length(grep("FINEMAP", p)) > 0 ) {
        data.table(
            file = p,
            data = sprintf("%s_%s", strsplit(p, "/")[[1]][3], strsplit(p, "/")[[1]][4]), 
            region = sub(".*sss/(.*)/.*", "\\1", p),
            chr = sub(".*sss/(.*):.*", "\\1", p),
            method = "FINEMAP-sss")
    }

}

# res_i = res[1]
# res_i = res[12]
# get_credset <- function(res_i) {
#     f <- res_i$file
#     method <- res_i$method

#     if (method == "BFMAP_FS") {
#         dat <- fread(f)
#         # dat[, .N, by = signal]
#         credset <- dat[, SNPname, by = signal]

#     } else if (method == "FINEMAP_sss") {
#         credp <- sub(".snp", ".cred1", f)
#         cred <- fread(credp, skip = 5)        
#         # cred[, sum(prob1)]
#         credset <- cred$cred1
#     }

#     return(credset)
# }

# overlap_perc_BFMAP <- function(credset_ori, credset) {
#     signals <- unique(credset$signal)
#     lapply(signals, function(s) {
#         x <- credset_ori[signal == s, SNPname]
#         y <- credset[signal == s, SNPname]
#         sum(y %in% x) / length(y)
#     })

# }

# overlap_perc_FINEMAP <- function(credset_ori, credset) {
#         sum(credset %in% credset_ori) / length(credset)
# }

process_file <- function(res_i) {
    f <- res_i$file
    dat <- fread(f)
    md <- res_i$method
    set <- res[method == md & chr == res_i$chr]

    if (md == "BFMAP-FS") {
        # credset_ori <- get_credset(set[QTL_window == "Original"])
        # credset <- get_credset(res_i)
        # op_perc <- overlap_perc_BFMAP(credset_ori, credset)

        leadsnp <- dat[order(-log_sBF, -normedProb)][
            1, .(leadSNP = SNPname, pos = Pos, log10bf = log_sBF, prob = normedProb)]

    } else if (md == "FINEMAP-sss") {
        # credset_ori <- get_credset(set[QTL_window == "Original"])
        # credset <- get_credset(res_i)
        # op_perc <- overlap_perc_FINEMAP(credset_ori, credset)

        leadsnp <- dat[order(-log10bf, -prob)][
            1, .(leadSNP = rsid, pos = position, log10bf = log10bf, prob = prob)]
    }

    return(leadsnp)
}



#

# testphenos <- c("BW32", "EW30", "EW70", "EN1", "EN13")

res <- rbindlist(lapply(filesp, extract_info))
res[, QTL_window := c("+1Mb", "Original", "-1Mb"), by = .(data, chr, method)]

res <- cbind(res, 
    rbindlist(lapply(1:nrow(res), function(i) {
    process_file(res[i,])
    }))
)
res[, dist_bp := pos - pos[QTL_window == "Original"], by = .(data, chr, method)]

res <- res[, -c("file", "region", "pos")][
    order(data, chr, method)
    ]

dataorder <- paste0(c("EN1", "EN13", "EW30", "EW70", "BW32"), "_all")
res$data <- factor(res$data, levels = dataorder)
setorder(res, data, chr)
write.table(res, outxt, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
