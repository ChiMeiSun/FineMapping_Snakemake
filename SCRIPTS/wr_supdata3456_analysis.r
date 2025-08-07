# args <- c("RESULTS/supdata/supdata.xlsx", 
# "RESULTS/BFMAPrecalc/BW32/geneprobs_BW32_all.txt")

# args <- c("RESULTS/supdata/supdata.xlsx", 
# "RESULTS/BFMAPrecalc/BW32/geneprobsori_BW32_all.txt")

# args <- c("RESULTS/supdata/supdata.xlsx", 
# "RESULTS/BFMAPrecalc/EN2/credsets_EN2_all.txt")

args <- commandArgs(TRUE)
args

library(data.table)
library(writexl)

files <- args[2:length(args)]
sheet_names <- tools::file_path_sans_ext(basename(files))
if (length(grep("enrichcat_", sheet_names)) > 0) {mark <- "ec"}
if (length(grep("geneprobs", sheet_names)) > 0) {mark <- "gp"}
if (length(grep("geneprobsori", sheet_names)) > 0) {mark <- "gpori"}
if (length(grep("credsets_", sheet_names)) > 0) {mark <- "cs"}

sheet_names <- sub("credsets_", "", sheet_names)
sheet_names <- sub("enrichcat_", "", sheet_names)
sheet_names <- sub("geneprobsori_", "", sheet_names)
sheet_names <- sub("geneprobs_", "", sheet_names)
dat <- setNames(lapply(files, fread), sheet_names)

# write_xlsx(dat, path = args[1])


### for analysis
if (mark == "gp" | mark == "gpori"){
    cc <- lapply(seq_along(dat), function(i){
            tmp <- dat[[i]]
            try({setnames(tmp, "summed_prob_re", "summed_prob")}, silent = TRUE)
            tmp <- tmp[summed_prob > 0.1 & gene_id != ""]
            tmp[gene_name == "", gene_name := gene_id]
            counts <- tmp[, .N, by = .(gene_name, chr)]
            counts[, dataset := names(dat)[i]]
            return(counts)
    })

    metacc <- rbindlist(cc, fill = TRUE)

    metacc$dataset <- factor(metacc$dataset, levels = unique(metacc$dataset))
    metacc_dc <- dcast(metacc, gene_name ~ dataset, value.var = "N", fill = 0)
    chrtab <- metacc[,unique(gene_name), chr]
    metacc_dc <- metacc_dc[chrtab, on = .(gene_name = V1)]
    if(mark == "gpori"){
        write.table(metacc_dc, "RESULTS/supdata/supdata3_c0.1.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
    }
    if(mark == "gp"){
        write.table(metacc_dc, "RESULTS/supdata/supdata6_c0.1.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
    }

}

if (mark == "cs"){
    cc <- lapply(seq_along(dat), function(i){
            tmp <- dat[[i]]
            tmp <- tmp[cat_group == "Consequence",]
            counts <- tmp[, .N, by = Chr]
            counts[, tot := sum(counts$N)]
            counts[, dataset := names(dat)[i]]
            return(counts)
    })

    metacc <- rbindlist(cc)
    metacc$dataset <- factor(metacc$dataset, levels = unique(metacc$dataset))
    write.table(metacc, "RESULTS/supdata/supdata5_cc.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

    cc <- lapply(seq_along(dat), function(i){
            tmp <- dat[[i]]
            tmp <- tmp[renormedProb > 0.9]
            counts <- tmp[, .N, by = Chr]
            counts[, dataset := names(dat)[i]]
            return(counts)
    })

    metacc <- rbindlist(cc)
    metacc$dataset <- factor(metacc$dataset, levels = unique(metacc$dataset))

    write.table(metacc, "RESULTS/supdata/supdata5_c0.9.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
}