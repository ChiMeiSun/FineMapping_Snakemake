# args <- c("RESULTS/gcta/mlma", "RESULTS/gcta/mlma_impute", 
# "RESULTS/PCA/array/IMAGE_QC.prune.in", "RESULTS/PCA/impute/bcftools_lifted_imputed_gt_QC.prune.in",
# "RESULTS/gcta/lambdaGC.txt")

args <- commandArgs(TRUE)
args

library(data.table)

pathar <- args[1]
pathimp <- args[2]
prunear <- args[3]
pruneimp <- args[4]
outxt <- args[5]

sufxmlma <- "mlma"

gen <- c('all','offspring','FounderBC1','BC1','FounderBC2','BC2','IC')
phenos <- c('EN1','EN2','EN3','EN4','EN5','EN6','EN7','EN8','EN10','EN11','EN12','EN13','BW32','EW30','EW40','EW50','EW70')
# phenos <- c('BW32')
gen <- factor(gen, levels = gen)
phenos <- factor(phenos, levels = phenos)

pruneinar <- fread(prunear, header = FALSE)
pruneinimp <- fread(pruneimp, header = FALSE)

cj <- CJ(gen = gen, pheno = phenos)
dat_arr <- copy(cj)[, file := sprintf("%s/%s/%s_%s.%s", pathar, gen, pheno, gen, sufxmlma)]
dat_imp <- copy(cj)[, file := sprintf("%s/%s/%s_%s.%s", pathimp, gen, pheno, gen, sufxmlma)]

dat_arr[, lambdaGC_arr := {
    p <- tryCatch({
        dat <- fread(file)
        dat[SNP %in% pruneinar$V1, p]
    }, error = function(e) return(NA)
    )
    median(qchisq(1-p, 1), na.rm = TRUE) / qchisq(0.5, 1)

}, by = file][, file := NULL]

dat_imp[, lambdaGC_imp := {
    p <- tryCatch({
        dat <- fread(file)
        dat[SNP %in% pruneinimp$V1, p]
    }, error = function(e) return(NA)
    )
    median(qchisq(1-p, 1), na.rm = TRUE) / qchisq(0.5, 1)
    
}, by = file][, file := NULL]

comb <- merge(dat_arr, dat_imp, by = c("gen", "pheno"), all = TRUE)

write.table(comb, outxt, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
