# install.packages("devtools")
# library(devtools)
# devtools::install_github("jiang18/gemrich")
# install.packages("ggnewscale")

# args = c("RESULTS/BFMAP/BW32/credsets_BW32_all.txt",      
# "RESULTS/BFMAP/BW32/vep_BW32_all.vcf",               
# "DATA/geno/ensembl_genes_GRCg6a.txt",                
# "RESULTS/BFMAPrecalc/BW32/credsets_BW32_all.txt",    
# "RESULTS/BFMAPrecalc/BW32/geneprobs_BW32_all.txt",   
# "RESULTS/BFMAPrecalc/BW32/enrichcat_BW32_all.txt",   
# "RESULTS/BFMAPrecalc/BW32/geneprobsori_BW32_all.txt",
# "PLOTS/BFMAP/BW32/recalc_BW32_all.pdf")
args <- commandArgs(TRUE)
args

library(gemrich)
library(data.table)
library(ggplot2)
library(patchwork) 
library(RColorBrewer)
library(ggnewscale)
library(stringr)

# ?estimate_category_enrichment
# ?bootstrap_category_enrichment
# ?renormalize_prob_by_enrichment
# ?calc_gene_posterior_prob
# ?map_snp_annotation
# ?calc_category_coverage
# ?calc_snp_category_prop





##### read data #####
credset <- fread(args[1])

vep <- fread(args[2], skip = "#Uploaded_variation")
colnames(vep)[1] <- "SNPname"
# check
# snps = credset[(duplicated(credset$SNPname)), SNPname]
# credset[SNPname %in% snps]

genes = fread(args[3])
colnames(genes) <- c("chr","start","end","gene_id","gene_id_ver","gene_name","strand","GOterm_acc","GOterm_name")
unigenes <- genes[, .SD[1], by = gene_id] # keep only unique variants

tragen <- sub(".*credsets_(.*)\\.txt", "\\1", args[1])
regions <- readLines(args[1], n = 1)
regions <- unlist(strsplit(sub("#", "", regions), ","))

svep <- vep[, .(SNPname, Consequence, Extra)]
svep[, IMPACT := str_extract(Extra, "(?<=IMPACT=)[^;]+")]
svep[, BIOTYPE := str_extract(Extra, "(?<=BIOTYPE=)[^;]+")]
svep <- svep[,-c("Extra")]

path = dirname(args[1])

cats <- c("Consequence", "IMPACT","BIOTYPE")
# VEP impact ranking: lower number = higher impact
impact_priority <- c(
  "splice_region_variant" = 1,
  "splice_polypyrimidine_tract_variant" = 2,
  "missense_variant" = 3,
  "synonymous_variant" = 4,
  "5_prime_UTR_variant" = 5,
  "3_prime_UTR_variant" = 6,
  "non_coding_transcript_exon_variant" = 7,
  "non_coding_transcript_variant" = 8,
  "intron_variant" = 9,
  "upstream_gene_variant" = 10,
  "downstream_gene_variant" = 11,
  "intergenic_variant" = 12
)

# plot & write table
renormed_bfmap <- list()
gene_probs <- list()
gene_probs_ori <- list()
mle <- list()
idx = 1

# parameters
min_snps = 10
min_p = 0.01
pvalue_threshold = 5e-5


make_mle_fail <- function(unisnpannot, cat, reason = 0) {
  msg_mle <<- paste0("Cat_group: ",cat,", Recalc PPC failed on estimate_category_enrichment") # assign to global variable
  warning("Make estimate_category_enrichment failed due to bad data")
  print(table(unisnpannot$multi_cat))

  msg_rea <- ifelse(reason == 1, paste0("num_var < ",min_snps),
              ifelse(reason == 2, paste0("type_cat < ",2), 
              ifelse(reason == 3, paste0("MLE error"), paste0(""))))

  res <- list()
  unicat <- unique(unisnpannot$multi_cat)
  if (length(unicat) == 0) unicat <- NA
  res$prob_mle <- data.table(category = unicat, q = NA, p_MLE = NA, Hessian_SE = NA)
  res$enrichment_mle <- data.table(category = unicat, enrichment = 1, SE = NA, p = NA, mle_status = "Failed", reason = msg_rea)
  return(res)
}

pdf(args[length(args)], width = 15, height = 9)

##### calculation by each chromosome #####
for (i in 1:length(regions)){
  reg = regions[i]
  c = as.numeric(unlist(strsplit(reg, ":"))[1])
  st = as.numeric(unlist(strsplit(reg, "[:-]"))[2])
  ed = as.numeric(unlist(strsplit(reg, "[:-]"))[3])
  
  print(paste0("Running QTL ",reg))

  file = paste0(path,"/VC/",reg,"_",tragen,".varcomp.csv")
  h2 = fread(file, nrows=2)[2,2]
  
  qtl_cs <- credset[Chr == c & Pos >= st & Pos <= ed,]
  kept_signals <- qtl_cs[SNPindex == 0 & Pval <= pvalue_threshold, signal]
  filter_cs <- qtl_cs[signal %in% kept_signals, ]

  # Calculate gene-level probabilities original
  gene_probs_c_ori <- calc_gene_posterior_prob(
    bfmap = qtl_cs,
    gene_annot = unigenes,
    extension = 1000, # extend gene boundaries by 1kb (total 2kb)
    pvalue_threshold = 1 # relaxed p-val to have a look of results before renormalization
  )
  gene_probs_c_ori <- gene_probs_c_ori[summed_prob >= 0.01,][order(-summed_prob)]
  gene_probs_c_ori[GOterm_name == "", GOterm_name := "/"]
  gene_probs_ori[[i]] <- gene_probs_c_ori

  gene_probs_c_ori[, mid := (start + end) / 2 ]
  gene_probs_c_ori[, signal := as.factor(signal)]
  randn <- rep(1:6, length = nrow(gene_probs_c_ori))
  gene_probs_c_ori$randn <- as.factor(randn)
  print(gene_probs_c_ori)
  
  ##### calculation with different categories #####
  for (cat in cats){ 
    print(paste0("Running QTL ",reg,", cat_group: ",cat))

    if (nrow(filter_cs) > 0) {
      snpannot <- svep[filter_cs[, .(Chr, SNPname)], on = .(SNPname = SNPname)]
        print("duplicated SNPs from vep")
        print(snpannot[duplicated(SNPname) | duplicated(SNPname, fromLast = TRUE)])
      setnames(snpannot, cat, "multi_cat")
      snpannot <- snpannot[,.(SNPname, multi_cat)]

      snpannot <- snpannot[, .(SNPname = rep(SNPname, sapply(strsplit(multi_cat, ","), length)),
                    multi_cat = unlist(strsplit(multi_cat, ",")))]
      
      if (cat == "Consequence"){
        ##### remove duplicated SNPs based on category ranking #####
        snpannot[, rank := impact_priority[multi_cat]]
        snpannot[is.na(rank), rank := 99]
        unisnpannot <- snpannot[order(rank), .SD[1], by = SNPname] # keep only unique variants
        # :SD -> subset. [1] means to keep only the first row of the subsetted group
        unisnpannot <- unisnpannot[!is.na(multi_cat)]
      } else{
        print(sum(duplicated(snpannot$SNPname)))
        unisnpannot <- snpannot
        unisnpannot <- unisnpannot[!is.na(multi_cat)]
      }

    } else {
      unisnpannot <- data.table(multi_cat = NA)
    }
    
    saveunisnpannot <- copy(unisnpannot)
    bfmap <- filter_cs[SNPname %in% unisnpannot$SNPname,]
    # if lead snp does not have annotation, and we still want MLE to run
    idr <- bfmap[, if (all(SNPindex != 0)) .I[which.min(Pval)], by = signal]$V1
    bfmap[idr, SNPindex := 0]

    
    ##### MLE category enrichment
    if (nrow(unisnpannot) < min_snps) {
      mle_result_c <- make_mle_fail(unisnpannot, cat, 1)
    } else if (length(unique(unisnpannot$multi_cat)) < 2) {
      mle_result_c <- make_mle_fail(unisnpannot, cat, 2)
    } else {
      # if cat prop too small, merge as "remaining"
      tab <- table(unisnpannot$multi_cat)
      cat_prop_c <- as.data.table(tab / nrow(unisnpannot))
      colnames(cat_prop_c) = c("category", "prop")
      sparse_cats <- cat_prop_c[prop < min_p, category]
      if (length(sparse_cats) == 1) {sparse_cats = NULL}
      unisnpannot[, multi_cat := ifelse(multi_cat %in% sparse_cats, "remaining", multi_cat)]
      
      tab <- table(unisnpannot$multi_cat)
      cat_prop_c <- as.data.table(tab / nrow(unisnpannot))
      colnames(cat_prop_c) = c("category", "prop")

      # Estimate enrichments using MLE
      # dont use return() unless inside a function
      mle_result_c <- tryCatch({
        res <- estimate_category_enrichment(
          bfmap = bfmap,
          snpinfo = unisnpannot,
          cat_prop = cat_prop_c,
          annot = "multi_cat",
          pvalue_threshold)
        
        msg_mle <<- paste0("Cat_group: ",cat)
        res$enrichment_mle[, mle_status := NA]
        res$enrichment_mle[, reason := NA]
        res
      },
      error = function(e) {
        warning(e$message)
        make_mle_fail(unisnpannot, cat, 3)
      })

    }

    mle_c <- data.table(mle_result_c$enrichment_mle,
                        mle_result_c$prob_mle[,3:4],
                        convergence = mle_result_c$convergence)[order(-enrichment)]
    mle_c[, data := tragen]
    mle_c[, QTL := reg]
    mle_c[, cat_group := cat]


    # Renormalize probabilities using enrichment estimates
    renormed_bfmap_c <- tryCatch({
      res <- renormalize_prob_by_enrichment(
              bfmap = bfmap,
              snpinfo = unisnpannot, 
              mle_result_c$enrichment_mle,
              annot = "multi_cat")
      res
      },
      error = function(e){
        cols <- c(colnames(qtl_cs), "category", "renormedProb","SNPindex2")
        res <- data.table(matrix(NA, nrow = 0, ncol = length(cols)))
        setnames(res, cols)
        res
    })
    renormed_bfmap_c <- renormed_bfmap_c[, -c("SNPindex2")]
    renormed_bfmap_c[, cat_group := cat]
    renormed_bfmap_c[, diff_prob := renormedProb - normedProb]
    renormed_bfmap_c <- renormed_bfmap_c[order(signal, -renormedProb)]

    # Calculate gene-level probabilities renormalized
    tmp <- renormed_bfmap_c[,c("signal", "SNPindex", "Chr", "Pos", "Pval", "renormedProb")]
    setnames(tmp,"renormedProb","normedProb")
    gene_probs_c_re <- calc_gene_posterior_prob(
      bfmap = tmp,
      gene_annot = unigenes,
      extension = 1000, # extend gene boundaries by 1kb
      pvalue_threshold
    )
    gene_probs_c_re[, mid := (start + end) / 2 ]
    gene_probs_c_re[, signal := as.factor(signal)]
    tmpid <- match(gene_probs_c_re$gene_name, gene_probs_c_ori$gene_name)
    gene_probs_c_re$randn <- as.factor(gene_probs_c_ori[tmpid, randn])
    gene_probs_c_re[GOterm_name == "", GOterm_name := "/"]

    if (nrow(gene_probs_c_re) > 0) {
      selcols <- setdiff(colnames(gene_probs_c_re), c("summed_prob","randn"))
      gene_probs_c <- merge(gene_probs_c_ori, gene_probs_c_re, by = selcols, all = TRUE, suffixes = c("_ori", "_re"))
      gene_probs_c[, randn_ori := NULL]
      gene_probs_c[, randn_re := NULL]

      gene_probs_c[, cat_group := cat]
      gene_probs_c[, diff_prob := summed_prob_re - summed_prob_ori]
      gene_probs_c <- gene_probs_c[order(signal, -summed_prob_re)]
    } else {
      cols <- c(setdiff(colnames(gene_probs_c_re), c("summed_prob", "randn")), 
                c("summed_prob_ori", "summed_prob_re", "cap_group", "diff_prob"))
      gene_probs_c <- data.table(matrix(NA, nrow = 0, ncol = length(cols)))
      colnames(gene_probs_c) <- cols
    }

    
    # plot
    theme = theme(text = element_text(size = 20), plot.margin = unit(c(0.5,0.5,1,1), "cm"), plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14))
    st <- as.numeric(unlist(strsplit(reg, "[:-]"))[2])
    ed <- as.numeric(unlist(strsplit(reg, "-"))[2])
    jitwid = (ed - st)/50
    jitwidzoom <- (max(gene_probs_c_ori$end) - min(gene_probs_c_ori$start)) / 50
    ptrhd = paste0("5e-5")

    renormed_bfmap_c[, category := as.factor(category)]
    renormed_bfmap_c[, signal := as.factor(signal)]

    colors <- scale_color_manual(values = c(brewer.pal(9, "Set1"), brewer.pal(12, "Paired")))
    
    # ppc vs recalc ppc
    ggppc <- ggplot(renormed_bfmap_c) +
              geom_point(aes(x = normedProb, y = renormedProb, color = category), size = 3) +
              labs(title = paste0(tragen,", QTL-",reg,", PPC vs recalc_PPC"), 
                    subtitle = paste0("filtered p-val: ",ptrhd,", ",msg_mle),
                    x = "ori_PPC", y = "recalc_PPC") +
              theme +
              colors

    # plot mle result
    ggmle = ggplot(mle_c) +
            geom_point(aes(x = p_MLE, y = enrichment, color = category), size = 3, alpha = 0.7) +
            labs(title = paste0(tragen,", QTL-",reg,", MLE of category"), 
                subtitle = paste0("filtered p-val: ",ptrhd,", ",msg_mle),
                x = "prob of inc causal SNP", y = "Enrichment", color = "category") +
            theme +
            colors

    # plot ppc per variant
    subgenes <- genes[chr == c & start >= st & end <= ed]
    subgenes <- subgenes[, .SD[1], by = gene_id] # remove duplicates
    subgenes[, mid := (start + end) / 2]
    randn <- rep(1:6, length = nrow(subgenes))
    subgenes$randn <- as.factor(randn)
    
    qtl_cs[, signal:= as.factor(signal)]
    ggop = ggplot(qtl_cs) +
            geom_point(aes(x = Pos, y = normedProb, color = signal), size = 2, alpha = 0.7) +
            labs(title = paste0(tragen,", QTL-",reg,", ori_PPC (no filter)"), 
                subtitle = paste0("h2: ",h2),
                x = "position", y = "ori_PPC") +
            ylim(-0.5,1) + xlim(st,ed) +
            colors +
            theme 
    ggrp = ggplot(renormed_bfmap_c) +
            geom_point(aes(x = Pos, y = renormedProb, color = signal), size = 2, alpha = 0.7) +
            labs(title = paste0(tragen,", QTL-",reg,", recalc_PPC"), 
                subtitle = paste0("filtered p-val: ",ptrhd,", ",msg_mle),
                x = "position", y = "recalc_PPC") +
            colors +

            new_scale_color() +
            geom_segment(data = subgenes, aes(x = start, y = -0.1, xend = end, yend = -0.1, color = randn),
                    arrow = arrow(length = unit(0.02, "npc")), inherit.aes=FALSE) +
            geom_text(data = subgenes, aes(x = mid, y = -0.3, label = gene_name, color = randn),
                size = 3, angle = 90, position = position_dodge(width=jitwid)) + 
            scale_color_manual(values = c(brewer.pal(9, "Set1"), brewer.pal(12, "Paired")), guide = "none") +
            ylim(-0.5,1) + xlim(st,ed) +
            theme 


    # plot ppc per gene
    ggeneori = ggplot(gene_probs_c_ori) +
            geom_point(aes(x = mid, y = summed_prob, color = signal), size = 2, alpha = 0.7) +
            geom_segment(aes(x = start, xend = end, y = -0.1, yend = -0.1, color = GOterm_name),
                    arrow = arrow(length = unit(0.05, "npc")), inherit.aes=FALSE, show.legend = FALSE) +
            labs(title = paste0(tragen,", QTL-",reg,", ori summed_prob of genes(zoom in)"), 
                subtitle = paste0("A brief look: no p-val filter & prob >= 0.01"),
                x = "Genes position", y = "sum_prob", color = "signal/GOterm") +
            geom_text(aes(x = mid, y = -0.3, label = gene_name, color = GOterm_name),
                size = 3, angle = 90, position = position_dodge(width=jitwidzoom)) + 
            colors +
            theme +
            ylim(-0.5,1)
    ggenere = ggplot(gene_probs_c_re) +
            geom_point(aes(x = mid, y = summed_prob, color = signal), size = 2, alpha = 0.7) +
            geom_segment(aes(x = start, xend = end, y = -0.1, yend = -0.1, color = GOterm_name),
                    arrow = arrow(length = unit(0.05, "npc")), inherit.aes=FALSE, show.legend = FALSE) +
            labs(title = paste0(tragen,", QTL-",reg,", recalc summed_prob of genes(zoom in)"), 
                subtitle = paste0("filtered p-val: ",ptrhd,", ",msg_mle),
                x = "Genes position", y = "sum_prob", color = "signal/GOterm") +
            geom_text(aes(x = mid, y = -0.3, label = gene_name, color = GOterm_name),
                size = 3, angle = 90, position = position_dodge(width=jitwidzoom)) + 
            colors +
            theme +
            ylim(-0.5,1)


    # save table
    gene_probs_c <- gene_probs_c[summed_prob_re >= 0.01,]
    if( !is.na(mle_c$mle_status[1])) { # if mle failed
      renormed_bfmap_c[, mle_status := "Failed"]
      gene_probs_c[, mle_status := "Failed"]
    } else{
      renormed_bfmap_c[, mle_status := NA]
      gene_probs_c[, mle_status := NA]
    }

    if (nrow(renormed_bfmap_c) > 0) {
      remainsnps <- renormed_bfmap_c[category == "remaining", SNPname]
      newcat <- paste0("remaining:", saveunisnpannot[renormed_bfmap_c[SNPname %in% remainsnps], on = .(SNPname), multi_cat])
      renormed_bfmap_c[SNPname %in% remainsnps, category := newcat]
    }

    mle[[idx]] <- mle_c
    renormed_bfmap[[idx]] <- renormed_bfmap_c
    gene_probs[[idx]] <- gene_probs_c
    idx = idx + 1

    # pdf("test.pdf", width = 15, height = 9)
    print(ggppc/ggmle)
    print(ggop/ggrp)
    print(ggeneori/ggenere)
  }
}
dev.off()

renormed_bfmap <- rbindlist(renormed_bfmap, fill = TRUE)
gene_probs <- rbindlist(gene_probs, fill = TRUE)
gene_probs_ori <- rbindlist(gene_probs_ori, fill = TRUE)
mle <- rbindlist(mle, fill = TRUE)
setcolorder(mle, c("data", "QTL", "cat_group", "category", "mle_status", "reason"))

write.table(renormed_bfmap, args[4], quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
write.table(gene_probs, args[5], quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
write.table(mle, args[6], quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t", na="")
write.table(gene_probs_ori, args[7], quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
