# args <- c("RESULTS/BFMAP/EW70/QTLs_EW70_all.txt", "RESULTS/gcta/mlma_impute/all/EW70_all.mlma",
# "INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC.fam", 
# "INTERMEDIATE/IDlist/IDlist_all.txt", "DATA/pheno/EW70.txt",
# "RESULTS/FINEMAP/EW70/all")

args <- commandArgs(TRUE)
args

#
library(data.table)
#
qtlsp <- args[1]
mlmap <- args[2]
famp <- args[3]
idp <- args[4]
phep <- args[5]
outpath <- args[6]
#

qtls <- fread(qtlsp)
mlma <- fread(mlmap)
fam <- fread(famp)
bfile <- gsub(".fam", "", famp)
idlist <- fread(idp, header = FALSE)
phe <- fread(phep)

aniids <- data.table(intersect(intersect(fam$V2, idlist$V2), phe$V2))
nsam <- nrow(aniids)

masterLD <- sprintf("%s/masterLD", outpath)
masterFM_sss <- sprintf("%s/sss/masterFM", outpath)
masterFM_cond <- sprintf("%s/cond/masterFM", outpath)

linesLD <- c("z;bcor;incl;ld;bgen;bgi;sample;n_samples")
linesFM_sss <- c("z;ld;snp;config;cred;log;n_samples")
linesFM_cond <- c("z;ld;snp;config;cred;log;n_samples")


for (i in seq_len(nrow(qtls))) {
    # test minor alteration of regions
    # tmp <- qtls[i]
    # tqtls <- data.table(Chr = tmp$Chr,
    #             st = c(tmp$exdst, tmp$exdst - 1e6, tmp$exdst + 1e6),
    #             ed = c(tmp$exded, tmp$exded + 1e6, tmp$exded - 1e6)
    # )
    tqtls <- qtls
    
    for (j in seq_len(nrow(tqtls))) {
    
        chr <- tqtls$Chr[j]
        st <- tqtls$st[j]
        ed <- tqtls$ed[j]
        region <- sprintf("%s:%s-%s", chr, format(st, scientific = FALSE), format(ed, scientific = FALSE))    
        
        namei <- sprintf("%s/%s/data", outpath, region)
        sssi <- sprintf("%s/sss/%s/data", outpath, region)
        condi <- sprintf("%s/cond/%s/data", outpath, region)

        lineldi <- paste(c(
                    paste0(namei, c(".z", ".bcor", ".incl", ".ld")), 
                    paste0(bfile, c(".bgen", ".bgen.bgi", ".sample")), 
                    nsam), collapse = ";"
                )
        linefmi_sss <- paste(c(
                paste0(namei, c(".z", ".ld")), 
                paste0(sssi, c(".snp", ".config", ".cred", ".log")), 
                nsam), collapse = ";"
            )
        linefmi_cond <- paste(c(
                paste0(namei, c(".z", ".ld")), 
                paste0(condi, c(".snp", ".config", ".cred", ".log")), 
                nsam), collapse = ";"
            )
        linesLD <- c(linesLD, lineldi)
        linesFM_sss <- c(linesFM_sss, linefmi_sss)
        linesFM_cond <- c(linesFM_cond, linefmi_cond)

        sub <- mlma[Chr == chr & bp >= st & bp <= ed,
            .(rsid = SNP, chromosome = Chr, position = bp, allele1 = A1, allele2 = A2, maf = Freq, beta = b, se)]

        # rm mkdir 
        fp <- sprintf("%s/%s", outpath, region)
        system(sprintf("rm -rf %s", fp))
        system(sprintf("mkdir %s", fp))

        fp <- sprintf("%s/sss/%s", outpath, region)
        system(sprintf("rm -rf %s", fp))
        system(sprintf("mkdir %s", fp))

        fp <- sprintf("%s/cond/%s", outpath, region)
        system(sprintf("rm -rf %s", fp))
        system(sprintf("mkdir %s", fp))

        write.table(sub, sprintf("%s.z", namei), sep = " ", quote = FALSE, col.names = TRUE, row.names = FALSE)
        write.table(aniids, sprintf("%s.incl", namei), sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
    }
}
writeLines(linesLD, masterLD)
writeLines(linesFM_sss, masterFM_sss)
writeLines(linesFM_cond, masterFM_cond)





