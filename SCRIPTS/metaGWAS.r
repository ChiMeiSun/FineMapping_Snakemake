# args <- c("RESULTS/gcta/mlma_impute","RESULTS/metaGWAS/ENs_all.txt","all")

args <- commandArgs(TRUE)
args

library(data.table)

phenos <- paste0("EN", c(1:8,10:13))
gen <- args[3]

pathimp <- args[1]


fgwa <- list()
file <- sprintf("%s/%s/%s_%s.mlma", pathimp, gen, phenos[1], gen)
save <- tryCatch(fread(file), error = function(e) return(NULL))
save <- save[,-c("p")]

for (p in phenos) {
    print(p)
    file <- sprintf("%s/%s/%s_%s.mlma", pathimp, gen, p, gen)
    dt <- tryCatch(fread(file, select = c("SNP","b","se")), error = function(e) return(NULL))
    if (!is.null(dt)) {
        dt[, t := b / se]
        fgwa[[p]] <- dt[,.(SNP, t)]
    }
}

snps <- fgwa[[1]]$SNP
tmat <- do.call(cbind, lapply(fgwa, function(x) x[match(snps, x$SNP), t]))
rownames(tmat) <- snps
invV <- solve(cor(tmat, use = "pairwise.complete.obs"))

chisq_vals <- rowSums((tmat %*% invV) * tmat) # get only diagonals
pvals <- pchisq(chisq_vals, df = ncol(tmat), lower.tail = FALSE)
if (length(pvals) != nrow(save)) warning("Num of SNPs do not match!")

metagwas <- cbind(save, p = pvals)
write.table(metagwas, args[2], row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

