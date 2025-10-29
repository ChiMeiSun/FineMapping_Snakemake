# args <- c("RESULTS/gcta/fastGWA","RESULTS/metaGWAS/ENlate_all.txt","all")

args <- commandArgs(TRUE)
args

library(data.table)

phenos <- paste0("EN", c(3:8,10:13))
gen <- args[3]
chrs = c(1:15,17:28,33)

pathimp <- args[1]


fgwa <- list()
file <- paste0(pathimp,"/",gen,"/",phenos[1],"_",gen,".fastGWA")
save <- tryCatch(fread(file), error = function(e) return(NULL))
colnames(save) = c("Chr", "SNP","bp","A1","A2","N","Freq","b","se","p")
save <- save[,-c("N","p")]

for (p in phenos) {
    print(p)
    file <- paste0(pathimp,"/",gen,"/",p,"_",gen,".fastGWA")
    dt <- tryCatch(fread(file, select = c("SNP","BETA","SE")), error = function(e) return(NULL))
    if (!is.null(dt)) {
        dt[, t := BETA/SE]
        fgwa[[p]] <- dt[,.(SNP, t)]
    }
}

snps <- fgwa[[1]]$SNP
tmat <- do.call(cbind, lapply(fgwa, function(x) x[match(snps, x$SNP), t]))
rownames(tmat) <- snps
invV <- solve(cor(tmat, use = "pairwise.complete.obs"))

chisq_vals <- rowSums((tmat %*% invV) * tmat) # get only diagonals
pvals <- pchisq(chisq_vals, df = ncol(tmat), lower.tail = FALSE)
length(pvals) == nrow(save)

metagwas <- cbind(save, p = pvals)
write.table(metagwas, args[2], row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

