
# args <- c("RESULTS/gcta/mlma_impute/all/BW32_all.mlma",                                   
# "DATA/geno/ensembl_genes_GRCg6a.txt",    
# "RESULTS/BFMAP/BW32/QTLs_BW32_all.txt",                                     
# "RESULTS/BFMAP/BW32/credsets_BW32_all.txt",        
# "RESULTS/BFMAP/BW32/h2_BW32_all.txt", 
# "PLOTS/BFMAP/BW32/BW32_all.pdf",    
# "INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC.bim",                 
# "RESULTS/BFMAP/BW32/FS/1:167395369-174490482_BW32_all.csv", 
# "RESULTS/BFMAP/BW32/FS/1:168395369-173490482_BW32_all.csv", 
# "RESULTS/BFMAP/BW32/FS/1:169395369-172490482_BW32_all.csv" #...    
# )


args <- commandArgs(TRUE)
args

library(data.table)
library(ggplot2)
library(patchwork) 

#
mlmap <- args[1]
genesp <- args[2]
qtlsp <- args[3]
outcred <- args[4]
outh2 <- args[5]
outplot <- args[6]
bimp <- args[7]
filesp <- args[8:length(args)]
#


mlma = fread(mlmap, header = TRUE)
dim(mlma)
bonf <- 0.05/nrow(mlma)


genes = fread(genesp, header = TRUE, na.string="")
colnames(genes) <- c("chr","start","end","gene_id","gene_id_ver","gene_name","strand","GOterm_acc","GOterm_name")
# unigenes <- genes[!is.na(chr)][, .SD[1], by = gene_id] # keep only unique variants

qtls <- fread(qtlsp)
qtls[, name := sprintf("%s:%s-%s", Chr, exdst, exded)]
# get the middle one if multiple QTLs
trueqtls <- qtls[, .SD[which.min(diff)], by = Chr]

fs <- lapply(filesp, function(p) {
    fread(p, header = TRUE)
})
names(fs) <- sub(".*FS/([^_]+)_.*", "\\1", filesp)
summary(fs)

credset <- rbindlist(fs[names(fs) %in% trueqtls$name], fill = TRUE)
header <- paste0("#", paste(names(fs)[names(fs) %in% trueqtls$name], collapse = ","))
writeLines(header, outcred)
write.table(credset, outcred, append = TRUE, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


tragen <- sub(".*QTLs_(.*).txt$", "\\1", qtlsp)

###

## save h2
sum_h2 <- list()
###

path <- dirname(qtlsp)

# plot
pdf(outplot, width = 15, height = 9)
break_ld <- c(0,0.2,0.4,0.6,0.8,1)
col_break <- c("royalblue","skyblue","seagreen","darkorange","red")
theme <- theme(text = element_text(size = 20), plot.margin = unit(c(0.5,0.5,1,1), "cm"))

for (i in 1:length(fs)){
    dat <- fs[[i]]
    dat[, signal := as.factor(signal)]
    region <- names(fs)[i]
    cc <- as.numeric(unlist(strsplit(region,"[:\\-]"))[1])
    st <- as.numeric(unlist(strsplit(region,"[:\\-]"))[2])
    ed <- as.numeric(unlist(strsplit(region,"[:\\-]"))[3])

    file <- sprintf("%s/VC/%s_%s.varcomp.csv", path, region, tragen)
    h2 <- fread(file, nrows = 2)[2,2][[1]]
    bim <- bimp
    cmd <- sprintf(
        "awk -v chrom=%s -v start=%s -v end=%s '$1 == chrom && $4 >= start && $4 <= end {print $2}' %s | wc -l",
        shQuote(cc),
        st,
        ed,
        shQuote(bim)
    )    
    nqtl <- system(cmd, intern = TRUE)
    sum_h2[[i]] <- list(tragen = tragen, region = region, h2 = h2, var_cs = length(unique(dat$SNPname)), var_qtl = nqtl)
    
    subgenes = genes[chr==cc & start>=st & end<=ed,]
    set.seed(123)
    if (length(unique(subgenes$gene_name)) <2){
    tab = data.table(gene_name = unique(subgenes$gene_name),
                    randn = abs(runif(1, 0, 1)),
                    tt = 2)
    } else {
            tab = data.table(gene_name = unique(subgenes$gene_name),
                            randn = abs(runif(unique(subgenes$gene_name), 0, 1)),
                            tt = rep(1:5, 1000)[1:length(unique(subgenes$gene_name))])
    }
    subgenes = tab[unique(subgenes, by="gene_name"), on="gene_name"]
    subgenes[is.na(gene_name), tt := 2] 
    subgenes[, tt := as.factor(tt)]
    subgenes[,mid := (start+end)/2]
    table(subgenes$tt)

    jitwid = (ed - st)/50
    maxbf = max(dat$log_sBF)
    minp = min(dat$Pval)
    g_bf = ggplot(dat) +
        geom_point(aes(x=Pos, y=log_sBF, color=(R)^2), size = 2, alpha = 0.5) +
        scale_y_continuous(name="log_sBF", limits=c(-2,maxbf+3)) +
        scale_color_gradientn(breaks = break_ld, colours = col_break, labels=break_ld, limits=c(0,1)) +
        labs(title = paste0("BFMAP, ",region,"_",tragen, ", h2 = ",round(h2,3)), 
            subtitle = paste0("Max log_sBF = ",maxbf,", min P-val = ",minp), x=paste0("chr",dat[1,Chr])) +
        theme + xlim(st,ed)
    
    g_prob = ggplot(dat) +
        geom_point(aes(x=Pos, y=normedProb, color=(R)^2), size = 2, alpha = 0.5) +
        scale_color_gradientn(breaks = break_ld, colours = col_break, labels=break_ld, limits=c(0,1)) +
        labs(title = paste0("BFMAP, ",region,"_",tragen, ", h2 = ",round(h2,3)), x=paste0("chr",dat[1,Chr])) +
        theme + 
        geom_segment(data = subgenes, aes(x=start, y=-0.1, xend=end, yend=-0.1, color=randn),
                arrow = arrow(length = unit(0.02, "npc")), arrow.fill = subgenes$randn, inherit.aes=FALSE) +
        geom_text(data = subgenes, aes(x = mid, y = -0.3, label = gene_name, color=randn),
            size = 3, angle = 90, position = position_dodge(width=jitwid)) + ylim(-0.5,1) + xlim(st,ed)

    g_cred = ggplot(dat) + 
        geom_point(aes(x = Pos, y = normedProb, group = signal, color = signal), size = 3, alpha = 0.7) +
        labs(title = paste0("BFMAP, ",region,"_",tragen, ", h2 = ",round(h2,3)), x=paste0("chr",dat[1,Chr])) +
        theme +
        geom_segment(data = subgenes, aes(x=start, y=-0.1, xend=end, yend=-0.1, color=tt),
                arrow = arrow(length = unit(0.02, "npc")), arrow.fill = subgenes$tt, inherit.aes=FALSE) +
        geom_text(data = subgenes, aes(x = mid, y = -0.3, label = gene_name, color=tt),
                size = 3, angle = 90, position = position_dodge(width=jitwid)) +ylim(-0.5, max(dat$normedProb)+0.2) + xlim(st,ed)

    m = mlma[Chr==cc & bp>=st & bp<=ed,]
    m = m[dat[,c("SNPname","R")], on=.(SNP = SNPname), nomatch=NULL]
    nrow(m) == nrow(dat)
    maxlp = -log10(min(m$p))
    g_fastgwa = ggplot(m) +
        geom_point(aes(x = bp, y = -log10(p), color = (R)^2), size = 2, alpha = 0.5) +
        scale_y_continuous(name = "-log10(p)", limits = c(0,maxlp+3)) +
        scale_color_gradientn(breaks = break_ld, colours = col_break, labels = break_ld, limits = c(0,1)) +
        geom_hline(yintercept = -log10(bonf), color = "black", linetype = "dashed") +
        labs(title = paste0("mlma, ",tragen), x = paste0("chr",dat[1,Chr])) +
        theme + xlim(st,ed)

        

print(g_bf/g_cred)
print(g_fastgwa/g_prob)
}
dev.off()

sum_h2 <- rbindlist(sum_h2)
write.table(sum_h2, outh2, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

