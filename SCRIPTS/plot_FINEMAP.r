# args <- c("RESULTS/FINEMAP/EW70/all/sss/masterFM", "RESULTS/FINEMAP/EW70/all/cond/masterFM",
# "DATA/geno/ensembl_genes_GRCg6a.txt",
# "RESULTS/gcta/mlma_impute/all/EW70_all.mlma",
# "PLOTS/FINEMAP/EW70/all/data.pdf"
# )

args <- commandArgs(TRUE)
args

library(data.table)
library(ggplot2)
library(patchwork) 
#
sssp <- args[1]
condp <- args[2]
genep <- args[3]
mlmap <- args[4]
outplot <- args[5]
#
qtls <- fread(sssp, select = 1)
qtls[, region := sub(".*all/(.*)/data.z$", "\\1", z) ]
qtls[, chr := sub(":.*", "", region) ]
qtls[, st := sub(".*:(.*)-.*", "\\1", region) ]
qtls[, ed := sub(".*-(.*)", "\\1", region) ]

sssp <- dirname(sssp)
condp <- dirname(condp)

tragen <- sub(".*all/(.*)\\.mlma$", "\\1", mlmap) 

mlma <- fread(mlmap, header=TRUE)
dim(mlma)
bonf <- 0.05/nrow(mlma)

genes <- fread(genep, header=TRUE, na.string="")
colnames(genes) = c("CHR","START","END","ENSG","ENSGver","GeneName","strand","GOaccs","GOname")



# plot
pdf(outplot, width = 15, height = 9)
ggformat <- list(
        theme(text = element_text(size = 20), plot.margin = unit(c(0.5,0.5,1,1), "cm")),
        scale_color_manual(values = c("0" = "steelblue", "1" = "red", 
                "2" = "grey", "3" = "darkorange", "4" = "seagreen", "5" = "purple")),
        scale_size_manual(values = c("0" = 2, "1" = 3.5))
)

for (i in seq_len(nrow(qtls))) {
        chr <- qtls$chr[i]
        st <- as.numeric(qtls$st[i])
        ed <- as.numeric(qtls$ed[i])
        region <- sprintf("%s:%s-%s", chr, format(st, scientific = FALSE), format(ed, scientific = FALSE))    

        # gene
        subgenes <- genes[CHR==chr & START>=st & END<=ed,]
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


        # --sss
        pp_sss <- sprintf("%s/%s", sssp, region) 

        snp_sss <- fread(sprintf("%s/data.snp", pp_sss))
        # dim = GWAS summary statistics
        # prob: marginal Posterior Inclusion Probabilities (PIP) of SNP causality
        # log10bf: log10 Bayes factors, quantifies the evidence that the l th SNP is causal with log10 Bayes factors greater than 2 reporting considerable evidence
        
        # config_sss <- fread(sprintf("%s/data.config", pp_sss))
        # prob: the posterior probabilities that configurations are the causal configuration
        # log10bf: log10 Bayes factors, quantifies the evidence for a causal configuration over the null configuration
        # h2: heritability contribution of SNPs
        # dim = 50k (default num of config to be saved; config can contain >1 SNP)
        
        msg_sss <- readLines(sprintf("%s/data.cred1", pp_sss), n = 5)
        cred1_sss <- fread(sprintf("%s/data.cred1", pp_sss), skip = 5)
        sum(cred1_sss$prob1)

        meta_sss <- merge(cred1_sss[, .(cred1)], 
                        snp_sss[, .(rsid, prob, log10bf)], by.x = "cred1", by.y = "rsid", all.x = TRUE)
        meta_sss[, color := 0]
        meta_sss[prob == max(prob) & log10bf == max(log10bf), color := 1]
        
        
        # --cond
        pp_cond <- sprintf("%s/%s", condp, region) 

        # snp_cond <- fread(sprintf("%s/data.snp", pp_cond))
        # dim = GWAS summary statistics

        msg_cond <- fread(sprintf("%s/data.config", pp_cond))
        condsnp <- msg_cond[rank == 2, config]
        # dim = 2 lines

        readLines(sprintf("%s/data.cred", pp_cond), n = 4)
        cred_cond <- fread(sprintf("%s/data.cred", pp_cond), skip = 4)
        sum(cred_cond$prob1)
        # dim = credset
        cred_cond[, color := 0]
        cred_cond[cred1 == condsnp, color := 1]

        # merge sss and cond
        colnames(meta_sss) <- paste0(colnames(meta_sss), "_sss")
        colnames(cred_cond) <- paste0(colnames(cred_cond), "_cond")
        meta <- merge(meta_sss, cred_cond, by.x = "cred1_sss", by.y = "cred1_cond", all = TRUE)
        meta <- merge(meta, mlma, by.x = "cred1_sss", by.y = "SNP", all.x = TRUE)

        ########
        jitwid = (ed - st)/50
        maxbf = max(meta$log10bf, na.rm = TRUE)
        minp = min(meta$p, na.rm = TRUE)

   g_bfs <- ggplot(meta) +
        geom_point(aes(x = bp, y = log10bf_sss, color = factor(color_sss), size = factor(color_sss)), alpha = 0.5) +
        ggformat +
        scale_y_continuous(name = "log_sBF", limits = c(-2, maxbf + 3)) +
        labs(title = paste0("FINEMAP_sss, ",region,"_",tragen), 
                subtitle = paste0("Max log_sBF = ",maxbf,", min P-val = ",minp), 
                x = paste0("chr",meta[1,Chr])) +
        xlim(st,ed)

   g_probs <- ggplot(meta) +
        geom_point(aes(x = bp, y = prob_sss, color = factor(color_sss), size = factor(color_sss)), alpha = 0.5) +
        ggformat +
        labs(title = paste0("FINEMAP_sss, ",region,"_",tragen), 
                subtitle = paste0(msg_sss[1],";",msg_sss[2]),
                x = paste0("chr",meta[1,Chr])) +
        geom_segment(data = subgenes, aes(x=START, y=-0.1, xend=END, yend=-0.1, color=tt),
                arrow = arrow(length = unit(0.02, "npc")), arrow.fill = subgenes$tt, inherit.aes=FALSE) +
        geom_text(data = subgenes, aes(x = mid, y = -0.3, label = GeneName, color=tt),
                size = 3, angle = 90, position = position_dodge(width=jitwid), inherit.aes=FALSE) + 
                ylim(-0.5,1) + xlim(st,ed)


   g_probc <- ggplot(meta) +
        geom_point(aes(x = bp, y = prob1_cond, color = factor(color_cond), size = factor(color_cond)), alpha = 0.5) +
        ggformat +
        labs(title = paste0("FINEMAP_cond, ",region,"_",tragen), x=paste0("chr",meta[1,Chr])) +
        geom_segment(data = subgenes, aes(x=START, y=-0.1, xend=END, yend=-0.1, color=tt),
                arrow = arrow(length = unit(0.02, "npc")), arrow.fill = subgenes$tt, inherit.aes=FALSE) +
        geom_text(data = subgenes, aes(x = mid, y = -0.3, label = GeneName, color=tt),
                size = 3, angle = 90, position = position_dodge(width=jitwid), inherit.aes=FALSE) + 
                ylim(-0.5,1) + xlim(st,ed)
                
        maxlp <- -log10(min(meta$p))
    g_mlma <- ggplot(meta) +
        geom_point(aes(x = bp, y = -log10(p), color = factor(color_sss), size = factor(color_sss)), alpha = 0.5) +
        ggformat +
        scale_y_continuous(name = "-log10(p)", limits = c(0, maxlp + 3)) +
        geom_hline(yintercept = -log10(bonf), color = "black", linetype = "dashed") +
        labs(title = paste0("mlma, ",tragen), x = paste0("chr", meta[1,Chr])) +
        xlim(st,ed)

        # pdf("test.pdf", width = 15, height = 9)
        print(g_bfs/g_probs)
        print(g_probc/g_mlma)
}

dev.off()
