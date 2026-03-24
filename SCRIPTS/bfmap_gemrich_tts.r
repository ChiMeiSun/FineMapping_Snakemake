# install.packages("devtools")
# library(devtools)
# devtools::install_github("jiang18/gemrich")
# install.packages("ggnewscale")


# args <- c("RESULTS/BFMAP/EN1/credsets_EN1_all.txt",      
# "RESULTS/BFMAP/EN1/vep_EN1_all.vcf",               
# "DATA/geno/ensembl_genes_GRCg6a.txt",                
# "RESULTS/BFMAPrecalc/ENs/credsets_ENs_all.txt",    
# "RESULTS/BFMAPrecalc/ENs/geneprobs_ENs_all.txt",   
# "RESULTS/BFMAPrecalc/ENs/enrichcat_ENs_all.txt",   
# "RESULTS/BFMAPrecalc/ENs/geneprobsori_ENs_all.txt",
# "PLOTS/BFMAP/ENs/recalc_ENs_all.pdf", 
# "ENs")

# args <- c("RESULTS/BFMAP/BW32/credsets_BW32_all.txt",      
# "RESULTS/BFMAP/BW32/vep_BW32_all.vcf",               
# "DATA/geno/ensembl_genes_GRCg6a.txt",                
# "RESULTS/BFMAPrecalc/BWEWs/credsets_BWEWs_all.txt",    
# "RESULTS/BFMAPrecalc/BWEWs/geneprobs_BWEWs_all.txt",   
# "RESULTS/BFMAPrecalc/BWEWs/enrichcat_BWEWs_all.txt",   
# "RESULTS/BFMAPrecalc/BWEWs/geneprobsori_BWEWs_all.txt",
# "PLOTS/BFMAP/BWEWs/recalc_BWEWs_all.pdf", 
# "BWEWs")

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


credp <- args[1]
vepp <- args[2]
genep <- args[3]
outcred <- args[4]
outgp <- args[5]
outerh <- args[6]
outgpori <- args[7]
outplot <- args[8]
tts <- args[9]


##### read data #####
genes = fread(genep)
colnames(genes) <- c("chr","start","end","gene_id","gene_id_ver","gene_name","strand","GOterm_acc","GOterm_name")
unigenes <- genes[, .SD[1], by = gene_id] # keep only unique variants


cats <- c("Consequence", "IMPACT","BIOTYPE")
# VEP impact ranking: lower number = higher impact
impact_priority <- c(
  "missense_variant" = 1,
  "splice_donor_5th_base_variant" = 2,
  "splice_region_variant" = 3,
  "splice_polypyrimidine_tract_variant" = 4,
  "synonymous_variant" = 5,
  "5_prime_UTR_variant" = 6,
  "3_prime_UTR_variant" = 7,
  "non_coding_transcript_exon_variant" = 8,
  "non_coding_transcript_variant" = 9,
  "intron_variant" = 10,
  "upstream_gene_variant" = 11,
  "downstream_gene_variant" = 12,
  "intergenic_variant" = 13
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


#
pdf(outplot, width = 15, height = 9)

# combine traits
regions <- readLines(credp, n = 1)
regions <- unlist(strsplit(sub("#", "", regions), ","))

path <- dirname(dirname(credp))

if (tts == "ENs") {
  tt <- c('EN1','EN2','EN3','EN4','EN5','EN6','EN7','EN8','EN10','EN11','EN12','EN13')

} else if (tts == "BWEWs") {
  tt <- c('BW32','EW30','EW40','EW50','EW70')

} else {
  stop("STOP: use 'ENs' to combine EN; use 'BWEWs' to combine BW and EWs; other names need to be defined!")
}

##### calculation by each chromosome #####
for (i in 1:length(regions)){
  reg = regions[i]
  c = as.numeric(unlist(strsplit(reg, ":"))[1])
  st = as.numeric(unlist(strsplit(reg, "[:-]"))[2])
  ed = as.numeric(unlist(strsplit(reg, "[:-]"))[3])
  
  tragen <- sub(".*recalc_(.*).pdf$", "\\1", outplot)
  if (c == 1 | c == 27) tragen <- sub("BWEWs", "BW32", tragen)

  print(sprintf("Running QTL %s for %s", reg, tragen))
  
  # rm(filter_cs)
  filter_cs <- rbindlist(lapply(tt, function (p) {

    cred_i <- fread(sprintf("%s/%s/credsets_%s_all.txt", path, p, p))
    cred_i <- cred_i[Chr == c & Pos >= st & Pos <= ed,]
    kept_signals <- cred_i[SNPindex == 0 & Pval <= pvalue_threshold, signal]
    filter_cs_i <- cred_i[signal %in% kept_signals, ]
    filter_cs_i[, pheno := p]
  }))



  # rm(qtl_cs)
  qtl_cs <- rbindlist(lapply(tt, function (p) {

    cred_i <- fread(sprintf("%s/%s/credsets_%s_all.txt", path, p, p))
    cred_i <- cred_i[Chr == c & Pos >= st & Pos <= ed,]
    cred_i[, pheno := p]
  }))
  qtl_cs <- qtl_cs[order(-normedProb), .SD[1], by = SNPname]

  # rewrite signal for it to run
  grp <- qtl_cs[, .GRP, by = .(pheno, signal)][, GRP := GRP - 1]
  
  qtl_cs <- qtl_cs[grp, on = .(pheno, signal)][
    ,savesignal := signal][, signal := GRP][, GRP := NULL]

  check <- qtl_cs[, .N, by = .(signal, SNPindex)][N>1]
  if (nrow(check) > 0) warning("Multiple signals have identical signal index")



  # rm(vep)
  vep <- rbindlist(lapply(tt, function (p) {

    vep_i <- fread(sprintf("%s/%s/vep_%s_all.vcf", path, p, p), skip = "#Uploaded_variation")
    colnames(vep_i)[1] <- "SNPname"
    vep_i <- vep_i[, .(SNPname, Consequence, Extra)]
    vep_i[, IMPACT := str_extract(Extra, "(?<=IMPACT=)[^;]+")]
    vep_i[, BIOTYPE := str_extract(Extra, "(?<=BIOTYPE=)[^;]+")]
    vep_i <- vep_i[,-c("Extra")]
    vep_i[, pheno := p]
  }))


  h2 = NA
  
  # Calculate gene-level probabilities original
  gene_probs_c_ori <- calc_gene_posterior_prob(
    bfmap = qtl_cs,
    gene_annot = unigenes,
    extension = 1000, # extend gene boundaries by 1kb (total 2kb)
    pvalue_threshold = 1 # relaxed p-val to have a look of results before renormalization
  )

  # combine duplicated genes prob across different signals (signals modified to fit the function thus not informative)
  # cols_keep <- setdiff(colnames(gene_probs_c_ori), c("summed_prob", "lead_pvalue"))
  # gene_probs_c_ori <- gene_probs_c_ori[, .(
  #   summed_prob = sum(summed_prob, na.rm = TRUE),
  #   lead_pvalue = min(lead_pvalue, na.rm = TRUE)
  # ), by = cols_keep]

  gene_probs_c_ori <- gene_probs_c_ori[summed_prob >= 0.01,]
  
  tab <- qtl_cs[, .N, by = .(signal, pheno, savesignal)][,N := NULL]
  gene_probs_c_ori <- merge(gene_probs_c_ori, tab, by = "signal", all.x = TRUE)
  
  gene_probs_c_ori <- gene_probs_c_ori[, signal := savesignal][
    GOterm_name == "", GOterm_name := "/"
  ][
    order(-summed_prob, pheno)
  ][, savesignal := NULL]

  gene_probs_ori[[i]] <- gene_probs_c_ori
  gene_probs_c_ori[, mid := (start + end) / 2 ]

  
  ##### calculation with different categories #####
  # cat <- "Consequence"
  for (cat in cats){ 
    print(paste0("Running QTL ",reg,", cat_group: ",cat))

    if (nrow(filter_cs) > 0) {
      snpannot <- vep[filter_cs[, .(Chr, SNPname)], on = .(SNPname = SNPname)]
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

    # rewrite signal for it to run
    grp <- filter_cs[, .GRP, by = .(pheno, signal)][, GRP := GRP - 1]
    
    filter_cs_rewrite <- copy(filter_cs)
    filter_cs_rewrite <- filter_cs_rewrite[grp, on = .(pheno, signal)][
    ,savesignal := signal][, signal := GRP][, GRP := NULL]

    bfmap <- filter_cs_rewrite[SNPname %in% unisnpannot$SNPname,]
    check <- bfmap[, .N, by = .(signal, SNPindex)][N>1]
    if (nrow(check) > 0) warning("Multiple signals have identical signal index")

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
      # tab <- table(unisnpannot$multi_cat)
      # cat_prop_c <- as.data.table(tab / nrow(unisnpannot))
      # colnames(cat_prop_c) = c("category", "prop")
      # sparse_cats <- cat_prop_c[prop < min_p, category]
      # if (length(sparse_cats) == 1) sparse_cats = NULL
      # unisnpannot[, multi_cat := ifelse(multi_cat %in% sparse_cats, "remaining", multi_cat)]
      
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
        cols <- c(setdiff(colnames(qtl_cs), c("savesignal", "SNPindex2")), 
                  "category", "renormedProb","cat_group", "diff_prob", "mle_status")
        res <- data.table(matrix(NA, nrow = 0, ncol = length(cols)))
        setnames(res, cols)
        res
    })


    # Calculate gene-level probabilities renormalized
    tmp <- renormed_bfmap_c[,c("signal", "SNPindex", "Chr", "Pos", "Pval", "renormedProb", "pheno")]
    setnames(tmp,"renormedProb","normedProb")
    gene_probs_c_re <- calc_gene_posterior_prob(
      bfmap = tmp,
      gene_annot = unigenes,
      extension = 1000, # extend gene boundaries by 1kb
      pvalue_threshold
    )
    
    # 
    tab <- renormed_bfmap_c[, .N, by = .(signal, pheno, savesignal)][,N := NULL]
    gene_probs_c_re <- merge(gene_probs_c_re, tab, by = "signal", all.x = TRUE)
    gene_probs_c_re[, signal := savesignal][, savesignal := NULL]
    gene_probs_c_re[GOterm_name == "", GOterm_name := "/"][
      , mid := (start + end) / 2
    ]
    #
    renormed_bfmap_c[, signal := savesignal][
      , savesignal := NULL][, SNPindex2 := NULL]
    renormed_bfmap_c[, cat_group := cat]
    renormed_bfmap_c[, diff_prob := renormedProb - normedProb]
    renormed_bfmap_c <- renormed_bfmap_c[order(signal, -renormedProb)]



    # gene_probs_c_ori[, pheno := as.factor(pheno)]
    # gene_probs_c_re[, pheno := as.factor(pheno)]


    if (nrow(gene_probs_c_re) > 0) {
      selcols <- setdiff(colnames(gene_probs_c_re), c("summed_prob", "lead_pvalue"))
      gene_probs_c <- merge(gene_probs_c_ori, gene_probs_c_re, by = selcols, all = TRUE, suffixes = c("_ori", "_re"))
      gene_probs_c[
        , lead_pvalue := min(lead_pvalue_ori, lead_pvalue_re, na.rm = TRUE)
        ][
          , lead_pvalue_ori := NULL
          ][
          , lead_pvalue_re := NULL
          ]

      gene_probs_c[, cat_group := cat]
      gene_probs_c[, diff_prob := summed_prob_re - summed_prob_ori]
      gene_probs_c <- gene_probs_c[order(-summed_prob_re)]
    
    } else {
      cols <- c(setdiff(colnames(gene_probs_c_ori), c("summed_prob")), 
                c("summed_prob_ori", "summed_prob_re", "cat_group", "diff_prob"))
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
    ggmle <- ggplot(mle_c) +
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
    ggop <- ggplot(qtl_cs) +
            geom_point(aes(x = Pos, y = normedProb, color = signal), size = 2, alpha = 0.7) +
            labs(title = paste0(tragen,", QTL-",reg,", ori_PPC (no filter)"), 
                subtitle = paste0("h2: ",h2),
                x = "position", y = "ori_PPC") +
            ylim(-0.5,1) + xlim(st,ed) +
            colors +
            theme 
    ggrp <- ggplot(renormed_bfmap_c) +
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
    ggeneori <- ggplot(gene_probs_c_ori) +
            geom_point(aes(x = mid, y = summed_prob, color = pheno), size = 2, alpha = 0.7) +
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

    ggenere <- ggplot(gene_probs_c_re) +
            geom_point(aes(x = mid, y = summed_prob, color = pheno), size = 2, alpha = 0.7) +
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

write.table(renormed_bfmap, outcred, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
write.table(gene_probs, outgp, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
write.table(mle, outerh, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t", na="")
write.table(gene_probs_ori, outgpori, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
