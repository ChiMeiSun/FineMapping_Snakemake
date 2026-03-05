calc_snp_category_prop <- function(snplist, bed, category_list = NULL) {
  # Ensure the input is a data.table
  setDT(bed)
  setDT(snplist)

  if (ncol(bed) < 4) {
    stop("Input 'bed' should have the following first four columns: chr, start, end, category.")
  }
  if(!all( c("chr", "pos") %in% colnames(snplist) )) {
    stop("Input 'snplist' needs to have the following two columns: chr and pos.")
  }

  cat_list = unique(bed[[4]])

  # Handle category_list
  if (is.null(category_list)) {
    if( length(cat_list) == 1 ) {
      category_list = cat_list
    } else {
      stop(paste(c("'category_list' is required when there are multiple categories in 'bed'"), collapse=" "))
    }
  } else {
    category_list = unique(category_list)
    if(any(!(category_list %in% cat_list))) {
      stop(paste(c("Categories in 'category_list' should be available in 'bed'"), collapse=" "))
    }
  }

  cat("There are", nrow(snplist), "SNPs in 'snplist'.\n")
  flush.console()

  if (!is.numeric(snplist$pos)) {
    stop("Error: column 'pos' should be numeric.")
  }

  snplist <- snplist[order(chr, pos)]
  snplist[, `:=`(pos2 = pos, id = .I)]

  bed = bed[, 1:4]
  setnames(bed, c("chr", "start", "end", "category"))
  bed$chr <- gsub("^chr", "", bed$chr, ignore.case = TRUE)
  bed <- bed[order(chr, start)]
  bed <- bed[category %chin% category_list]

  snplist[, chr := as.character(chr)]
  bed[, chr := as.character(chr)]

  cat("Setting keys for 'snplist' and 'bed'...\n")
  flush.console()

  # Perform a join using data.table
  setkey(snplist, chr, pos, pos2)
  setkey(bed, chr, start, end)

  cat("Finding overlaps of 'snplist' and 'bed'...\n")
  flush.console()

  # Join SNPs to genomic features
  overlaps <- foverlaps(
    snplist,
    bed,
    by.x = c("chr", "pos", "pos2"),
    by.y = c("chr", "start", "end"),
    nomatch = 0
  )

  if (nrow(overlaps) == 0) {
    cat("No SNPs are within genomic features. No SNP info file generated.\n")
    return(NULL)
  }
  cat(length(unique(overlaps$id)), "SNPs are within genomic annotations.\n")
  flush.console()

  categories <- unique(overlaps$category)
  snp_info <- data.table(SNP = snplist$id)

  # Add categories with NA by default
  snp_info[, (categories) := NA_real_]

  # Create unique id-category pairs for efficiency
  unique_pairs <- unique(overlaps[, .(id, category)])
  # Set values using data.table syntax
  for (cat in categories) {
      matched_ids <- unique_pairs[category == cat, id]
      set(snp_info, i = which(snp_info$SNP %in% matched_ids), j = cat, value = 1)
  }

  if (length(categories) == 1){
    snp_info_merged = snp_info
    colnames(snp_info_merged)[2] = "multi_cat"
    snp_info_merged$multi_cat[is.na(snp_info_merged$multi_cat)] = "remaining"
  }else{
    # Create output data.table
    snp_info_merged = data.table(SNP = snp_info$SNP, multi_cat = NA_character_)  # Changed to character
    
    # Process all categories at once
    for (cat in category_list) {
      if(!(cat %in% categories)) next
      # Use set() for more efficient updates
      set(snp_info_merged, 
          i = which(snp_info[[cat]] == 1 & is.na(snp_info_merged$multi_cat)), 
          j = "multi_cat",
          value = cat)
    }
    snp_info_merged$multi_cat[is.na(snp_info_merged$multi_cat)] = "remaining"
  }

  result = table(snp_info_merged$multi_cat) / nrow(snp_info_merged)
  if (length(categories) == 1) names(result)[names(result) == "1"] = categories[1]
  result = as.data.table(result[c(category_list, "remaining")])
  colnames(result) = c("category", "prop")
  result[,1] = c(category_list, "remaining")
  result[is.na(prop), prop := 0]
  return(result[])
}


estimate_category_enrichment <- function(bfmap, snpinfo, cat_prop, annot = "multi_cat", pvalue_threshold = 5e-5) {
  # Validate inputs
  pvalue_threshold = as.numeric(pvalue_threshold)
  if(is.na(pvalue_threshold)) {
	stop("'pvalue_threshold' must be numeric.\n", call.=FALSE)
  }
  # Load data
  kept_signals <- bfmap$signal[bfmap$SNPindex == 0 & bfmap$Pval < pvalue_threshold]
  if(length(kept_signals) != length(unique(kept_signals))) {
    stop("Multiple signals have identical signal index.", call.=FALSE)
  }
  bfmap <- bfmap[bfmap$signal %in% kept_signals, ]
  setDT(bfmap)
  nloci = length(kept_signals)

  if (!(annot %in% colnames(snpinfo))) stop(paste(annot, "is not one of the column names in 'snpinfo'"))

  cat_prop[[1]] = factor(cat_prop[[1]], levels=cat_prop[[1]])
  cat_names = levels(cat_prop[[1]])

  snpinfo <- copy(as.data.table(snpinfo))
  snp_colname = colnames(snpinfo)[1]
  snpinfo = copy(snpinfo)[, .SD, .SDcols = c(snp_colname, annot)]
  snpinfo = snpinfo[snpinfo[[snp_colname]] %in% unique(bfmap$SNPname), ]
  snpinfo[is.na(snpinfo[[annot]]),2] = "remaining"
  if(! all(unique(snpinfo[[annot]]) %in% cat_names) ) {
    stop(paste("Some categories in 'snpinfo' are missing in 'cat_prop'.\n"), call.=FALSE)
  }
  snpinfo[[annot]] = factor(snpinfo[[annot]], levels=cat_names)
  cat_cnt = table(snpinfo[[annot]])
  
  levels(snpinfo[[annot]])= c(1:length(cat_names))
  var2cat = new.env(hash=TRUE)
  cats = apply(snpinfo, 1, function(x) var2cat[[x[1]]] = x[annot])
  bfmap[, cat_idx := as.numeric(unlist(sapply(SNPname, function(x) var2cat[[x]])))]

  cat("Number of unique model SNPs in each category:")
  print(cat_cnt)
  if(any(cat_cnt < 50)) {
    cat("\nWarning: Some categories involve too few signal SNPs, which may cause ML estimation problems.\n")
  }
  flush.console()

  if(any(cat_prop[[2]]<=0)) {
    stop(paste("Category proportions in cat_prop must be positive.\n"), call.=FALSE)
  }
  freq = cat_prop[[2]]
  if(sum(freq) != 1) {
    cat("\nSum of proportions in cat_prop is not equal to 1. Scaled to 1.\n")
    freq = freq / sum(freq)
  }

  cat(paste("\nCompleted reading all data files.\n",
      "The ML estimation uses ", nloci, " loci and ", nrow(bfmap)," variants.\n",sep=""))
  flush.console()

  logLL = function(par) {
    par <- c(par, 1-sum(par))   # sum(par) = 1
    if(!any(par<0)) {
      # li = rep(0, nloci)
      # for(i in c(1:nloci)) {
      #   d = bfmap[bfmap$signal==kept_signals[i],]
      #   pp = par[d$cat_idx]
      #   qp = freq[d$cat_idx]
      #   li[i] = sum(d[[pcol]]*pp/qp)
      # }
      ret = bfmap[, .(li = sum(normedProb * par[cat_idx] / freq[cat_idx])), by = signal][, sum(log(li))]
      return(ret)
    } else return(-100000)
  }

  par = freq[1:(length(freq)-1)]
  npars = length(par)
  result = NA

  eps = 1e-12
  if(npars == 1) {
    result <- optim(par, logLL, method="Brent",
            hessian=T,
            lower=rep(eps, npars),
            upper=rep(1-eps, npars),
            control=list(fnscale=-1))
  } else {
    # First try optim
    result <- try(optim(par, logLL, method="L-BFGS-B",
              hessian=T,
              lower=rep(eps, npars),
              upper=rep(1-eps, npars),
              control=list(fnscale=-1)), silent=TRUE)
    # If optim fails, try constrOptim
    if(inherits(result, "try-error") || result$convergence != 0 || any(result$par < 1e-3)) {
      cat("Initial optimization failed. Trying constrOptim...\n")
      flush.console()

      ui = rbind(diag(npars), rep(-1, npars))
      ci = c(rep(eps, npars), -1)
      result_constr <- constrOptim(theta=par, f=logLL, grad=NULL,
                  ui=ui, ci=ci,
                  control=list(fnscale=-1))
      # Add Hessian to result_constr before assigning to result
      result_constr$hessian <- try(optimHess(result_constr$par, logLL))
      result <- result_constr
      cat("constrOptim completed.\n")
      flush.console()
    }
  }
  cat("Completed MLE.\n")

  covar_mle = solve(-result$hessian)
  cat("Calculated SE with Hessian.\n")
  flush.console()

  # Profile likelihood method for SE
  profile_se <- function(mle, delta = seq(-0.1, 0.1, length=41)) {
    # Calculate profile likelihood around MLE
    profile_points <- mle + delta
    profile_points <- profile_points[profile_points>0 & profile_points<1]
    profile_values <- sapply(profile_points, logLL)

    # Maximum likelihood at MLE
    max_ll <- max(profile_values)

    # Find points where log-likelihood drops by 1.92
    # (approximately 95% CI region)
    cutoff <- max_ll - 1.92
    valid_points <- profile_points[profile_values >= cutoff]

    # SE estimate from the range
    se <- diff(range(valid_points))/3.84

    # Return full information for diagnostics
    return(list(
      se = se,
      profile_points = profile_points,
      profile_values = profile_values,
      cutoff = cutoff,
      ci = range(valid_points)
    ))
  }

  mle <- result$par
  prob_mle = copy(cat_prop)
  setDT(prob_mle)
  prob_mle[, paste0("V", 3:7) := NA_real_]
  prob_mle[[3]]=c(mle, 1-sum(mle))

  prob_cov_matrix <- rbind(
    cbind(covar_mle, -covar_mle %*% rep(1, nrow(covar_mle))),
    c(-covar_mle %*% rep(1, nrow(covar_mle)), sum(covar_mle))
  )
  prob_mle[[4]] = sqrt(diag(prob_cov_matrix))
  if(npars == 1) {
    profileLL_result <- profile_se(mle)
    cat("Completed profile likelihood (used only for binary snpinfoations).\n")
    prob_mle[[5]] = rep(profileLL_result$se, 2)
    prob_mle[1, (6:7) := as.list(profileLL_result$ci)]
    prob_mle[2, (6:7) := as.list(rev(1-profileLL_result$ci))]
  }
  colnames(prob_mle)[2:7] = c("q", "p_MLE", "Hessian_SE", "profile_SE", "profile_95lower", "profile_95upper")
  colnames(prob_cov_matrix) = levels(cat_prop[[1]])

  enrichment_mle = prob_mle[,1:4]
  colnames(enrichment_mle) =c("category","enrichment","SE","p")
  enrichment_mle[,2] = prob_mle$p_MLE / prob_mle$q
  enrichment_mle[,3] = prob_mle$Hessian_SE / prob_mle$q
  enrichment_mle[,4] = pnorm( (enrichment_mle[[2]]-1)/enrichment_mle[[3]], lower.tail = F)

  output <- list()
  output[["prob_mle"]] = prob_mle
  output[["prob_cov_matrix"]] = prob_cov_matrix
  output[["loglik"]] = result$value
  output[["convergence"]] = result$convergence
  output[["counts"]] = cat_cnt
  output[["enrichment_mle"]] = enrichment_mle
  return(output)
}


renormalize_prob_by_enrichment <- function(bfmap, snpinfo, enrichment_mle, annot = "multi_cat") {
  # Column checks
  required_bfmap <- c("signal", "SNPname", "normedProb")
  
  if (!all(required_bfmap %in% colnames(bfmap))) {
    stop("Missing bfmap columns: ", paste(setdiff(required_bfmap, colnames(bfmap)), collapse=", "))
  }
  if (!(annot %in% colnames(snpinfo))) stop(paste(annot, "is not one of the column names in 'snpinfo'"))

  bfmap <- copy(bfmap)
  setDT(bfmap)
  signals = unique(bfmap$signal)
  nloci = length(signals)

  enrichment_mle[[1]] = factor(enrichment_mle[[1]], levels=enrichment_mle[[1]])
  cat_names = levels(enrichment_mle[[1]])

  snp_colname = colnames(snpinfo)[1]
  snpinfo = snpinfo[,mget(c(snp_colname, annot))]
  snpinfo = snpinfo[snpinfo[[snp_colname]] %in% unique(bfmap$SNPname), ]
  snpinfo[is.na(snpinfo[[annot]]),2] = "remaining"
  if(! all(unique(snpinfo[[annot]]) %in% cat_names) ) {
    stop(paste("Some categories in 'snpinfo' are missing in 'enrichment_mle'.\n"), call.=FALSE)
  }
  snpinfo[[annot]] = factor(snpinfo[[annot]], levels=cat_names)
  
  levels(snpinfo[[annot]])= c(1:length(cat_names))
  var2cat = new.env(hash=TRUE)
  cats = apply(snpinfo, 1, function(x) var2cat[[x[1]]] = x[annot])
  bfmap[, cat_idx := as.numeric(unlist(sapply(SNPname, function(x) var2cat[[x]])))]

  if(any(enrichment_mle[[2]]<=0)) {
    stop(paste("Fold enrichments in enrichment_mle must be positive.\n"), call.=FALSE)
  }
  fold = enrichment_mle[[2]]

  cat(paste("\nCompleted reading all data files.\n",
      "Renormalization will cover ", nloci, " loci and ", nrow(bfmap)," variants.\n",sep=""))
  flush.console()
  
  bfmap[, category := cat_names[cat_idx]]
  bfmap[, renormedProb := {
    ff <- fold[cat_idx]
    normedProb * ff / sum(normedProb * ff) * sum(normedProb)
  }, by = signal]
  bfmap[, cat_idx := NULL]

  cat("Renormalization completed normally.\n")
  
  bfmap[order(-renormedProb), SNPindex2 := seq_len(.N)-1L, by=signal]
  setorder(bfmap, signal, SNPindex2)
  return(bfmap[])
}


calc_gene_posterior_prob <- function(bfmap, gene_annot, extension = 0, pvalue_threshold = 5e-5) {
  setDT(bfmap)
  setDT(gene_annot)
  
  bfmap[, Chr := as.character(Chr)]
  message("Starting calculations..."); flush.console()
  
  # Column checks
  required_bfmap <- c("signal", "SNPindex", "Chr", "Pos", "Pval", "normedProb")
  required_gene <- c("chr", "start", "end", "gene_id")
  
  if (!all(required_bfmap %in% colnames(bfmap))) {
    stop("Missing bfmap columns: ", paste(setdiff(required_bfmap, colnames(bfmap)), collapse=", "))
  }
  if (!all(required_gene %in% colnames(gene_annot))) {
    stop("Missing gene_annot columns: ", paste(setdiff(required_gene, colnames(gene_annot)), collapse=", "))
  }
  
  # Process reverse strand
  message("Processing gene coordinates..."); flush.console()
  reverse_idx <- gene_annot[, which(end < start)]
  if (length(reverse_idx) > 0) {
    message(sprintf("Found %d reverse strand genes", length(reverse_idx))); flush.console()
    gene_annot[reverse_idx, c("start", "end") := .(end, start)]
  }
  
  # Filter signals and extend boundaries
  message("Filtering signals..."); flush.console()
  sig_signals <- bfmap[SNPindex == 0 & Pval <= pvalue_threshold, unique(signal)]
  message(sprintf("Found %d significant signals", length(sig_signals))); flush.console()
  
  bfmap_filtered <- bfmap[signal %in% sig_signals]
  gene_annot[, `:=`(
    start_ext = start - extension,
    end_ext = end + extension
  )]
  
  # Calculate overlaps and sum probabilities
  message("Calculating overlaps..."); flush.console()
  results <- bfmap_filtered[gene_annot, on = .(Chr = chr), 
    allow.cartesian = TRUE, nomatch = 0][
      Pos >= start_ext & Pos <= end_ext,
      .(summed_prob = sum(normedProb)), 
      by = .(gene_id, signal)
    ]
  
  # Add gene information and lead p-values
  results <- gene_annot[results, on = "gene_id"]
  lead_pvals <- bfmap_filtered[SNPindex == 0, .(signal, lead_pvalue = Pval)]
  results <- lead_pvals[results, on = "signal"]
  
  message(sprintf("Found %d gene-signal pairs", nrow(results))); flush.console()

  setcolorder(results, c("gene_id", "chr", "start", "end", "signal", "summed_prob", "lead_pvalue", "start_ext", "end_ext"))
  setorder(results, -summed_prob, chr, start, signal)
  return(results[])
}

args = c("RESULTS/BFMAP/BW32/credsets_BW32_all.txt",
"RESULTS/BFMAP/BW32/vep_BW32_all.vcf",
"DATA/geno/ensembl_genes_GRCg6a.txt",
"RESULTS/BFMAPrecalc/BW32/credsets_BW32_all.txt",
"RESULTS/BFMAPrecalc/BW32/geneprobs_BW32_all.txt",
"RESULTS/BFMAPrecalc/BW32/enrichcat_BW32_all.txt",
"PLOTS/BFMAP/BW32/recalc_BW32_all.pdf")

args = c("RESULTS/BFMAP/EN2/credsets_EN2_all.txt",
"RESULTS/BFMAP/EN2/vep_EN2_all.vcf",
"DATA/geno/ensembl_genes_GRCg6a.txt",
"RESULTS/BFMAPrecalc/EN2/credsets_EN2_all.txt",
"RESULTS/BFMAPrecalc/EN2/geneprobs_EN2_all.txt",
"RESULTS/BFMAPrecalc/EN2/enrichcat_EN2_all.txt",
"PLOTS/BFMAP/EN2/recalc_EN2_all.pdf")

# args = c("../proj_IMAGE/test_bfmap/Outputs/bcf_impute/KGW/ForwardSelection_27:4000000-8000000_KGW.csv",
# "RESULTS/BFMAP/BW32/vep_BW32_all.vcf",
# "DATA/geno/ensembl_genes_GRCg6a.txt",
# "RESULTS/BFMAPrecalc/BW32/credsets_BW32_all.txt",
# "RESULTS/BFMAPrecalc/BW32/geneprobs_BW32_all.txt",
# "RESULTS/BFMAPrecalc/BW32/enrichcat_BW32_all.txt",
# "PLOTS/BFMAP/BW32/recalc_BW32_all.pdf")
# reg = "27:4000000-8000000"
# h2 = 0.402

# args <- c("RESULTS/gcta/fastGWA/IC/EW30_IC.fastGWA", "RESULTS/metaGWAS/ENlate_all.txt")
