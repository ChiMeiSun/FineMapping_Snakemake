
# args = c("RESULTS/gcta/fastGWA/all/EN2_all.fastGWA",                                   
# "DATA/geno/ensembl_genes_GRCg6a.txt",                                         
# "RESULTS/BFMAP/EN2/credsets_EN2_all.txt",                              
# "RESULTS/BFMAP/EN2/FS/1:26659374-29808199_EN2_all.csv",    
# "RESULTS/BFMAP/EN2/FS/1:40490453-43119109_EN2_all.csv"  ,
# "PLOTS/BFMAP/EN2/EN2_all.pdf")


args <- commandArgs(TRUE)
args

library(data.table)
library(ggplot2)
library(patchwork) 

fgwa = fread(args[1], header=TRUE)
dim(fgwa)
bonf = 0.05/nrow(fgwa)
# fgwa = fgwa[,-6]
# colnames(fgwa) = c("Chr", "SNP","bp","A1","A2","Freq","b","se","p")
# head(fgwa)

genes = fread(args[2], header=TRUE, na.string="")
colnames(genes) = c("CHR","START","END","ENSG","ENSGver","GeneName","strand","GOaccs","GOname")


ind_fs = grep("FS",args)
fs = list()
j=1
for (i in ind_fs){
    fs[[j]] = fread(args[i], header = TRUE)
    region = sub(".*FS/([^_]+)_.*", "\\1", args[i])
    names(fs)[j] = region
    j = j+1
}
summary(fs)

credset <- tryCatch({rbindlist(fs, fill = TRUE)}, error = function(e){warning("fuck")})
header <- paste0("#", paste(names(fs), collapse = ","))
writeLines(header, args[3])
write.table(credset, args[3], append=TRUE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

x = unlist(strsplit(args[1], "[/_]"))
trait = x[5]
gen = x[4]
tragen = paste0(trait,"_",gen)

###

## save h2
sum_h2 <- list()
###

path = dirname(args[3])

# plot
pdf(args[length(args)], width = 15, height = 9)
break_ld = c(0,0.2,0.4,0.6,0.8,1)
col_break = c("royalblue","skyblue","seagreen","darkorange","red")
theme = theme(text = element_text(size = 20), plot.margin = unit(c(0.5,0.5,1,1), "cm"))

for (i in 1:length(fs)){
    dat = fs[[i]]
    dat[,signal:=as.factor(signal)]
    region = names(fs)[i]
    cc = as.numeric(unlist(strsplit(region,"[:\\-]"))[1])
    st = as.numeric(unlist(strsplit(region,"[:\\-]"))[2])
    ed = as.numeric(unlist(strsplit(region,"[:\\-]"))[3])

    file = paste0(path,"/VC/",region,"_",tragen,".varcomp.csv")
    h2 = fread(file, nrows=2)[2,2][[1]]
    bim = "INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC.bim"
    cmd <- sprintf(
        "awk -v chrom=%s -v start=%s -v end=%s '$1 == chrom && $4 >= start && $4 <= end {print $2}' %s | wc -l",
        shQuote(cc),
        st,
        ed,
        shQuote(bim)
    )    
    nqtl = system(cmd, intern = TRUE)
    sum_h2[[i]] <- list(tragen = tragen, region = region, h2 = h2, var_cs = nrow(dat), var_qtl = nqtl)
    
    subgenes = genes[CHR==cc & START>=st & END<=ed,]
    set.seed(123)
    if (length(unique(subgenes$GeneName)) <2){
    tab = data.table(GeneName = unique(subgenes$GeneName),
                    randn = abs(runif(1, 0, 1)),
                    tt = 2)
    } else {
            tab = data.table(GeneName = unique(subgenes$GeneName),
                            randn = abs(runif(unique(subgenes$GeneName), 0, 1)),
                            tt = rep(1:5, 1000)[1:length(unique(subgenes$GeneName))])
    }
    subgenes = tab[unique(subgenes, by="GeneName"), on="GeneName"]
    subgenes[is.na(GeneName), tt := 2] 
    subgenes[, tt := as.factor(tt)]
    subgenes[,mid := (START+END)/2]
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
        geom_segment(data = subgenes, aes(x=START, y=-0.1, xend=END, yend=-0.1, color=randn),
                arrow = arrow(length = unit(0.02, "npc")), arrow.fill = subgenes$randn, inherit.aes=FALSE) +
        geom_text(data = subgenes, aes(x = mid, y = -0.3, label = GeneName, color=randn),
            size = 3, angle = 90, position = position_dodge(width=jitwid)) + ylim(-0.5,1) + xlim(st,ed)

    g_cred = ggplot(dat) + 
        geom_point(aes(x = Pos, y = normedProb, group = signal, color = signal), size = 3, alpha = 0.7) +
        labs(title = paste0("BFMAP, ",region,"_",tragen, ", h2 = ",round(h2,3)), x=paste0("chr",dat[1,Chr])) +
        theme +
        geom_segment(data = subgenes, aes(x=START, y=-0.1, xend=END, yend=-0.1, color=tt),
                arrow = arrow(length = unit(0.02, "npc")), arrow.fill = subgenes$tt, inherit.aes=FALSE) +
        geom_text(data = subgenes, aes(x = mid, y = -0.3, label = GeneName, color=tt),
                size = 3, angle = 90, position = position_dodge(width=jitwid)) +ylim(-0.5, max(dat$normedProb)+0.2) + xlim(st,ed)

    m = fgwa[CHR==cc & POS>=st & POS<=ed,]
    m = m[dat[,c("SNPname","R")], on=.(SNP = SNPname), nomatch=NULL]
    nrow(m) == nrow(dat)
    maxlp = -log10(min(m$P))
    g_fastgwa = ggplot(m) +
        geom_point(aes(x = POS, y = -log10(P), color = (R)^2), size = 2, alpha = 0.5) +
        scale_y_continuous(name = "-log10(p)", limits = c(0,maxlp+3)) +
        scale_color_gradientn(breaks = break_ld, colours = col_break, labels = break_ld, limits = c(0,1)) +
        geom_hline(yintercept = -log10(bonf), color = "black", linetype = "dashed") +
        labs(title = paste0("fastGWAS, ",tragen), x = paste0("chr",dat[1,Chr])) +
        theme + xlim(st,ed)

        

print(g_bf/g_cred)
print(g_fastgwa/g_prob)
}
dev.off()

sum_h2 <- rbindlist(sum_h2)
write.table(sum_h2, args[4], quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

