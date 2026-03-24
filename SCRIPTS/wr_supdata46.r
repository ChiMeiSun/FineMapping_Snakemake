# args <- c("RESULTS/supdata/supdata.xlsx", 
# "RESULTS/BFMAPrecalc/ENs/geneprobs_ENs_all.txt",
# "RESULTS/BFMAPrecalc/BWEWs/geneprobs_BWEWs_all.txt")

# args <- c("RESULTS/supdata/supdata.xlsx", 
# "RESULTS/BFMAPrecalc/BWEWs/credsets_BWEWs_all.txt")

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

if (mark == "gp") {
    dat <- lapply(dat, function(d) {
        d <- d[, -c("start_ext", "end_ext", "mid")]
        setorder(d, chr, gene_id, -summed_prob_re)
    })
}
write_xlsx(dat, path = args[1])
