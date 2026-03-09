args <- c("RESULTS/gcta/", "RESULTS/PCA/array/IMAGE_QC.prune.in", "RESULTS/PCA/impute/bcftools_lifted_imputed_gt_QC.prune.in")

library(data.table)

path <- args[1]
prunearp <- args[2]
pruneimp <- args[3]

sufxmlma <- "mlma"
sufxfastgwa <- "fastGWA"
sufxloco <- "loco.mlma"

gen <- c('all','offspring','FounderBC1','BC1','FounderBC2','BC2','IC')
phenos <- c('EN1','EN2','EN3','EN4','EN5','EN6','EN7','EN8','EN10','EN11','EN12','EN13','BW32','EW30','EW40','EW50','EW70')
phenos <- c('BW32')
gen <- factor(gen, levels = gen)
phenos <- factor(phenos, levels = phenos)

pruneinar <- fread(prunearp, header = FALSE)
pruneinimp <- fread(pruneimp, header = FALSE)

cj <- CJ(gen = gen, pheno = phenos)
dat_arml <- copy(cj)[, file := sprintf("%s/%s/%s/%s_%s.%s", path, sufxmlma, gen, pheno, gen, sufxmlma)]
dat_imfg <- copy(cj)[, file := sprintf("%s/%s/%s/%s_%s.%s", path, sufxfastgwa, gen, pheno, gen, sufxfastgwa)]

dat_arml[, lambdaGC := {
    dat <- fread(file)
    p <- dat[SNP %in% pruneinar$V1][[ncol(dat)]]
    # p <- dat[[ncol(dat)]]

    median(qchisq(1-p, 1), na.rm = TRUE) / qchisq(0.5, 1)

}, by = file]
mean(dat_arml$lambdaGC)

dat_imfg[, lambdaGC := {
    dat <- fread(file)
    p <- dat[SNP %in% pruneinimp$V1][[ncol(dat)]]
    p <- dat[[ncol(dat)]]

    median(qchisq(1-p, 1), na.rm = TRUE) / qchisq(0.5, 1)
    
}, by = file]
mean(dat_imfg$lambdaGC)
