# args <- c("RESULTS/FINEMAP/EW70/all/sss/masterFM", "RESULTS/FINEMAP/EW70/all/cond/masterFM",
# "DATA/geno/ensembl_genes_GRCg6a.txt",
# "RESULTS/gcta/mlma_impute/all/EW70_all.mlma",
# "RESULTS/FINEMAP/EW70/all/sss/credsets.txt",
# "RESULTS/FINEMAP/EW70/all/sss/geneprobs.txt",
# "PLOTS/FINEMAP/EW70/all/data.pdf")

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
outcs_sss <- args[5]
outgp_sss <- args[6]
outplot <- args[7]
#


qtls <- fread(sssp, select = 1)
qtls[, region := sub(".*all/(.*)/data.z$", "\\1", z) ]
qtls[, chr := sub(":.*", "", region) ]
qtls[, st := as.numeric(sub(".*:(.*)-.*", "\\1", region)) ]
qtls[, ed := sub(".*-(.*)", "\\1", region) ]

trueqtls <- qtls[, .SD[which(abs(st - median(st)) == min(abs(st - median(st))))[1]], by = chr]

sssp <- dirname(sssp)
condp <- dirname(condp)

tragen <- sub(".*all/(.*)\\.mlma$", "\\1", mlmap) 

mlma <- fread(mlmap, header=TRUE)
dim(mlma)
bonf <- 0.05/nrow(mlma)

genes <- fread(genep, header=TRUE, na.string="")
colnames(genes) <- c("chr","start","end","gene_id","gene_id_ver","gene_name","strand","GOterm_acc","GOterm_name")
unigenes <- genes[!is.na(chr)][, .SD[1], by = gene_id] # keep only unique variants

#
get_gene_prob <- function(dat, subgenes, extension = 1000) {
    subgenes[, `:=`(start_ext = start - extension, 
                        end_ext = end + extension)]
    # dat should contain: chr, position, prob, log10bf
    colnames(dat) <- c("chr", "position", "prob", "log10bf")
    dat[, chr := as.character(chr)]
    subgenes[, chr := as.character(chr)]
    
    res <- dat[subgenes, 
                   on = .(chr = chr, 
                   position >= start_ext,
                   position <= end_ext), 
                   nomatch = NULL]
    cols <- setdiff(colnames(subgenes), c("start_ext", "end_ext"))
    res <- res[, .(sum_prob = sum(prob, na.rm = TRUE),
                N = .N,
                maxlog10bf = max(log10bf, na.rm = TRUE)),
                by = cols]
    return(res)
}

#

# plot
pdf(outplot, width = 15, height = 9)
ggformat <- list(
        theme(text = element_text(size = 20), plot.margin = unit(c(0.5,0.5,1,1), "cm")),
        scale_color_manual(values = c("0" = "steelblue", "1" = "red", 
                "2" = "grey", "3" = "darkorange", "4" = "seagreen", "5" = "purple")),
        scale_size_manual(values = c("0" = 2, "1" = 3.5))
)

cs_sss <- data.table()
gp_sss <- data.table()
msg_sss <- c()

for (i in seq_len(nrow(qtls))) {
        chr <- qtls$chr[i]
        st <- as.numeric(qtls$st[i])
        ed <- as.numeric(qtls$ed[i])
        region <- sprintf("%s:%s-%s", chr, format(st, scientific = FALSE), format(ed, scientific = FALSE))    

        # gene
        subgenes <- unigenes[chr==chr & start>=st & end<=ed,]
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
        
        cred1_sss <- fread(sprintf("%s/data.cred1", pp_sss), skip = 5)
        sum(cred1_sss$prob1)

        meta_sss <- merge(cred1_sss[, .(cred1)], 
                        snp_sss[, .(rsid, prob, log10bf)], by.x = "cred1", by.y = "rsid", all.x = TRUE)
        meta_sss[, color := 0]
        meta_sss[prob == max(prob) & log10bf == max(log10bf), color := 1]
        
        ## save true qtl output table
        if (region %in% trueqtls$region) {

                msg_sss <- c(msg_sss,
                        sprintf("#QTL: %s", region),
                        readLines(sprintf("%s/data.cred1", pp_sss), n = 5)
                )
                cred2p <- sprintf("%s/data.cred2", pp_sss)
                if (file.exists(cred2p)) msg_sss <- c(msg_sss, readLines(cred2p, n = 5))

                cs_sss <- rbind(cs_sss,
                        merge(cred1_sss[, .(SNPanme = cred1)], 
                                snp_sss[, .(SNPanme = rsid, Chr = chromosome, Pos = position, prob, log10bf)], 
                                by = "SNPanme", all.x = TRUE)
                )
                
                gp_sss <- rbind(gp_sss,
                        get_gene_prob(cs_sss[, .(Chr, Pos, prob, log10bf)], subgenes, extension = 1000)
                )
        }
        ##

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


        # # get gene prob
        # ## sss
        # dat <- copy(meta[, .(chr = Chr, position = bp, prob = prob_sss, log10bf = log10bf_sss)])
        # gp_sss <- get_gene_prob(dat, subgenes, extension = 1000)
        # ## cond
        # dat <- copy(meta[, .(chr = Chr, position = bp, prob = prob1_cond, log10bf = NA)])
        # gp_cond <- get_gene_prob(dat, subgenes, extension = 1000)

        ######## plot ########
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
        geom_segment(data = subgenes, aes(x=start, y=-0.1, xend=end, yend=-0.1, color=tt),
                arrow = arrow(length = unit(0.02, "npc")), arrow.fill = subgenes$tt, inherit.aes=FALSE) +
        geom_text(data = subgenes, aes(x = mid, y = -0.3, label = gene_name, color=tt),
                size = 3, angle = 90, position = position_dodge(width=jitwid), inherit.aes=FALSE) + 
                ylim(-0.5,1) + xlim(st,ed)


   g_probc <- ggplot(meta) +
        geom_point(aes(x = bp, y = prob1_cond, color = factor(color_cond), size = factor(color_cond)), alpha = 0.5) +
        ggformat +
        labs(title = paste0("FINEMAP_cond, ",region,"_",tragen), x=paste0("chr",meta[1,Chr])) +
        geom_segment(data = subgenes, aes(x=start, y=-0.1, xend=end, yend=-0.1, color=tt),
                arrow = arrow(length = unit(0.02, "npc")), arrow.fill = subgenes$tt, inherit.aes=FALSE) +
        geom_text(data = subgenes, aes(x = mid, y = -0.3, label = gene_name, color=tt),
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

writeLines(msg_sss, outcs_sss)
write.table(cs_sss, outcs_sss, append=TRUE, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
write.table(gp_sss, outgp_sss, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")



# For --sss (sum of single effects):
# File	What it contains	How to use it
# data.snp	Marginal inclusion probabilities (PIP) for each SNP across all causal configurations	Primary file for candidate variants - This is what you should use for identifying likely causal SNPs
# data.cred (or .cred1, .cred2)	95% credible sets - groups of SNPs that together have 95% probability of containing the causal variant	Use to define credible sets for each signal. The sum of probabilities in each .cred file = 0.95
# data.config	Top ranked configurations (combinations of SNPs) with their posterior probabilities	Shows which combinations of SNPs are most likely. Each row is a different model
# data.bf (or .bf1, .bf2)	Bayes factors for each independent signal	Used to compare evidence for 1 vs 2 vs 3 causal variants