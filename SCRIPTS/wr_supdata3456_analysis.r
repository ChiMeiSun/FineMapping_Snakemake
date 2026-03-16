
library(data.table)
options(scipen=999, width=2000)

phenos = c('EN1','EN2','EN3','EN4','EN5','EN6','EN7','EN8','EN10','EN11','EN12','EN13','BW32','EW30','EW40','EW50','EW70')
phenos = c('ENs','BWEW')
datype = c("enrichcat","geneprobs","geneprobsori","credsets")
datype = c("geneprobs")

for (i in seq_along(datype)) {
    dat <- datype[i]

    files <- sprintf("RESULTS/BFMAPrecalc/%s/%s_%s_all.txt",phenos,dat,phenos)


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

    # for a gene
    gene <- c("FAM124A","ENSGALG00000052822","ENSGALG00000053256")
    gene <- c("IGF2BP1","GIP","SNF8","UBE2Z","ENSGALG00000001525")
    gene <- c("NCAPG","LCORL","LDB2","IGF2BP1")
    gene <- c("ENSGALG00000049955","GARS")
    gene <- c("XIRP1","OXSR1","DLEC1","ENSGALG00000005884","CTDSPL")
    gene <- c("ENSGALG00000050406","ENSGALG00000049302")
    gene <- c("ENSGALG00000005934","ENSGALG00000006018","OXSR1") # EN1 post-enrich
    gene <- c("LDB2","ENSGALG00000032384","SLIT2","ENSGALG00000054285","TAPT1", "PROM1", "FBXL5", "CD38")

    probs <- data.table()
    probs <- lapply(seq_along(dat), function(i){
        tmp <- dat[[i]]
        if (mark == "gpori") {
            tmp <- tmp[gene_id %in% gene | gene_name %in% gene,
                    .(gene_id, p = summed_prob, gene_name)]
        } else if (mark == "gp" & nrow(tmp) > 0) {
            tmp <- tmp[gene_id %in% gene | gene_name %in% gene,
                    .(gene_id, p = summed_prob_re, gene_name, cat_group, mle_status,summed_prob_ori)]
        } else{
            tmp <- data.table()
        }
        tmp[, pheno := names(dat[i])]
        return(tmp)
    })
    probs <- rbindlist(probs, fill = TRUE)
    probs <- probs[!is.na(gene_id)]
    probs[, .SD[order(-p)], by = gene_id]
    probs[is.na(mle_status), .SD[order(-p)], by = gene_id]
    probs[, .SD[order(-p, cat_group)], by = .(gene_id,pheno)]
    probs[, .SD[order(-p)]]

}

## compare h2
path = "RESULTS/gcta/reml"
gen = "all"
tt = c('EN1','EN2','EN3','EN4','EN5','EN6','EN7','EN8','EN10','EN11','EN12','EN13','BW32','EW30','EW40','EW50','EW70')
datah2 = data.table(path = paste0(path,"/",tt,"/",gen,".hsq"), pheno = tt, h2 = NA)

for ( i in 1:nrow(datah2)){
    tmp = fread(datah2$path[i], fill=TRUE)
    datah2$h2[i] = tmp[Source == "V(G)/Vp",Variance]
}
datah2[, path := NULL]
qtlh2 <- fread("RESULTS/supdata/SupplementaryData2.txt")
qtlh2[, tragen := sub("_all","",tragen)]
x = qtlh2[datah2, on=.(tragen = pheno)]
cor(x[,.(h2, i.h2)])