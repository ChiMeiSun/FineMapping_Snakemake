# args <- c("RESULTS/supdata/supdata2.txt", 
# "RESULTS/BFMAP/BW32/h2_BW32_all.txt",
# "RESULTS/BFMAP/EW40/h2_EW40_all.txt",
# "RESULTS/FINEMAP/BW32/all/sss/credsets.txt",
# "RESULTS/FINEMAP/EN1/all/sss/credsets.txt"
# )

args <- commandArgs(TRUE)
args

library(data.table)
library(ggplot2)

#
outxt <- args[1]
filesp <- args[2:length(args)]
tragenorder <- paste0(c('EN1','EN2','EN3','EN4','EN5','EN6','EN7','EN8','EN10','EN11','EN12','EN13','BW32','EW30','EW40','EW50','EW70'), "_all")
#
library(data.table)

#  function
parse_finemap_summary <- function(dat) {
  qtl_indices <- grep("^#QTL: ", dat)
  
  dt_list <- lapply(qtl_indices, function(idx) {
    # Get region
    region <- gsub("^#QTL: ", "", dat[idx])
    
    # Find next lines
    next_lines <- dat[(idx+1):min(idx+11, length(dat))]
    
    # Extract pr1
    pr1_line <- grep("^# Post-Pr\\(# of causal SNPs is 1\\)", next_lines, value = TRUE)
    pr1 <- if (length(pr1_line) > 0) as.numeric(gsub(".*= (.*)", "\\1", pr1_line[1])) else NA
    
    # Extract pr2
    pr2_line <- grep("^# Post-Pr\\(# of causal SNPs is 2\\)", next_lines, value = TRUE)
    pr2 <- if (length(pr2_line) > 0) as.numeric(gsub(".*= (.*)", "\\1", pr2_line[1])) else NA
    
    # Extract mean_ld1
    mean_line <- grep("^#mean\\(\\|ld\\|\\)", next_lines, value = TRUE)
    if (length(mean_line) > 0) {
      values <- (unlist(strsplit(gsub("^#mean\\(\\|ld\\|\\) ", "", mean_line), " ")))
      mean_ld1 <- values[1]
      mean_ld2 <- if (length(values) >= 3) sprintf("%s,%s",values[3], values[5]) else NA
    } else {
      mean_ld1 <- mean_ld2 <- NA
    }
    
    data.table(
      region = region,
      pr1 = pr1,
      pr2 = pr2,
      mean_ld1 = mean_ld1,
      mean_ld2 = mean_ld2
    )
  })
  
  rbindlist(dt_list)
}
#


filesp_bfmap <- filesp[grep("BFMAP", filesp)]
filesp_finemap <- filesp[grep("FINEMAP", filesp)]

res_bfmap <- rbindlist(lapply(filesp_bfmap, function (p) {

    dat <- fread(p)
    dat[, chr := sub(":.*", "", region) ]
    dat[, st := as.numeric(sub(".*:(.*)-.*", "\\1", region)) ]

    dat <- dat[, .SD[which(abs(st - median(st)) == min(abs(st - median(st))))[1]], by = chr]
    dat[, st := NULL]
}))

res_finemap <- rbindlist(lapply(filesp_finemap, function (p) {

    tragen <- sprintf("%s_%s", strsplit(p, "/")[[1]][3], strsplit(p, "/")[[1]][4] )
    dat <- readLines(p, 50)
    dat1 <- parse_finemap_summary(dat)
    dat1[, Chr := as.numeric(sub(":.*", "", region)) ]

    dat <- fread(p)
    dat2 <- dat[, .(FINEMAP_var_cs1 = .N), by = Chr]

    out <- merge(dat2, dat1, by = "Chr")
    out[, tragen := tragen][, Chr := NULL]
}))

# merge, rename
allres <- merge(res_bfmap, res_finemap, by = c("tragen", "region") )
setcolorder(allres, c("tragen", "chr", "region", "var_qtl"))
allres$tragen <- factor(allres$tragen, levels = tragenorder)
setorder(allres, tragen, chr)

cols_bfmap <- c("h2", "var_cs")
cols_finemap <- c("pr1", "pr2", "mean_ld1", "mean_ld2")
setnames(allres, cols_bfmap, paste0("BFMAP_", cols_bfmap))
setnames(allres, cols_finemap, paste0("FINEMAP_", cols_finemap))


write.table(allres, outxt, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

png(sub(".txt", ".png", outxt), width = 5000, height = 3000, res = 300)
ggplot(allres, aes(x = tragen, y = BFMAP_h2)) + geom_boxplot()
dev.off()
