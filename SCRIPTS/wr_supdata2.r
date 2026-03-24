# args <- c("RESULTS/supdata/SupplementaryData2.xlsx", 
# "RESULTS/analysis/BFMAPvsFINEMAP.txt",
# "RESULTS/BFMAPrecalc/BWEWs/credsets_BWEWs_all.txt",
# "RESULTS/BFMAPrecalc/ENs/credsets_ENs_all.txt",
# "RESULTS/BFMAP/BW32/h2_BW32_all.txt",
# "RESULTS/BFMAP/EW40/h2_EW40_all.txt",
# "RESULTS/FINEMAP/BW32/all/sss/credsets.txt",
# "RESULTS/FINEMAP/EN1/all/sss/credsets.txt"
# )

args <- commandArgs(TRUE)
args

library(data.table)
library(writexl)
library(ggplot2)

#
outp <- args[1]
f1p <- args[2]
dat1p <- args[3]
dat2p <- args[4]
filesp <- args[5:length(args)]
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

    # Extract log10bf1,2
    bf_line <- grep("^#log10bf ", next_lines, value = TRUE)
    if (length(bf_line) > 0) {
      values <- (unlist(strsplit(gsub("^#log10bf ", "", bf_line), " ")))
      bf1 <- values[1]
      bf2 <- if (length(values) >= 3) sprintf("%s,%s",values[3], values[5]) else NA
    } else {
      bf1 <- bf2 <- NA
    }
   
        
    # Extract mean_ld1,2
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
      log10bf1 = bf1,
      log10bf2 = bf2,
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
    dat[, chr := as.numeric(sub(":.*", "", region)) ]
    dat[, st := as.numeric(sub(".*:(.*)-.*", "\\1", region)) ]

    dat <- dat[, .SD[which(abs(st - median(st)) == min(abs(st - median(st))))[1]], by = chr]
    dat[, st := NULL]
}))

dat1 <- fread(dat1p)
dat2 <- fread(dat2p)
ndat1 <- dat1[cat_group == "Consequence", .N, by = .(Chr, pheno)][, pheno := paste0(pheno,"_all")]
ndat2 <- dat2[cat_group == "Consequence", .N, by = .(Chr, pheno)][, pheno := paste0(pheno,"_all")]

res_bfmap[ndat1, var_cs_filtered := i.N, on = .(tragen = pheno, chr = Chr)]
res_bfmap[ndat2, var_cs_filtered := i.N, on = .(tragen = pheno, chr = Chr)]
res_bfmap[is.na(var_cs_filtered), var_cs_filtered := 0]

#
res_finemap <- rbindlist(lapply(filesp_finemap, function (p) {

    tragen <- sprintf("%s_%s", strsplit(p, "/")[[1]][3], strsplit(p, "/")[[1]][4] )
    dat <- readLines(p, 50)
    dat1 <- parse_finemap_summary(dat)
    dat1[, Chr := as.numeric(sub(":.*", "", region)) ]

    dat <- fread(p)
    dat2 <- dat[, .(
      FINEMAP_var = .N,
      FINEMAP_var_evi = sum(log10bf >= 2, na.rm = TRUE)
      ), by = Chr]

    out <- merge(dat2, dat1, by = "Chr")
    out[, tragen := tragen][, Chr := NULL]
}))

# merge, rename
allres <- merge(res_bfmap, res_finemap, by = c("tragen", "region") )
setcolorder(allres, c("tragen", "chr", "region", "var_qtl"))
allres$tragen <- factor(allres$tragen, levels = tragenorder)
setorder(allres, tragen, chr)

cols_bfmap <- c("h2", "var_cs", "var_cs_filtered")
cols_finemap <- c("pr1", "pr2", "mean_ld1", "mean_ld2")
setnames(allres, cols_bfmap, paste0("BFMAP_", cols_bfmap))
setnames(allres, cols_finemap, paste0("FINEMAP_", cols_finemap))

# write excel
out <- list(fread(f1p), allres)
out <- setNames(out, paste0("Result", 1 : length(out)))
write_xlsx(out, path = outp)

