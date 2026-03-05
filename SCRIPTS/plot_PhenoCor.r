# args <- c("DATA/pheno/Phenotypes_Final.csv",
#   "PLOTS/publish/PhenoCor.png",
# "RESULTS/gcta/bireml/Imputed", "RESULTS/gcta/reml")

args <- commandArgs(TRUE)
args

library(data.table)
library(corrplot)

pheno <- fread(args[1])
colnames(pheno)
tt = c('EN1','EN2','EN3','EN4','EN5','EN6','EN7','EN8','EN10','EN11','EN12','EN13','BW32','EW30','EW40','EW50','EW70')


#### get h2 ####
path = args[4]
gen = "all"
datah2 = data.table(path = paste0(path,"/",tt,"/",gen,".hsq"), pheno = tt, h2 = NA)

for ( i in 1:nrow(datah2)){
    tmp = fread(datah2$path[i], fill=TRUE)
    datah2$h2[i] = tmp[Source == "V(G)/Vp",Variance]
}


#### get genetic correlations ####
  path = args[3]

  datagc = data.table(file = system(paste0("ls ",path," | grep hsq"), intern=TRUE),
                      t1 = NA, t2 = NA, cor_g = NA, se_cor_g = NA, h2_t1 = NA, h2_t2 = NA
  )
  datagc[, t1:= substr(file,1,4)]
  datagc[, t2:= substr(file,5,20)]
  datagc[, t1:= sub("_","",t1)]
  datagc[, t2:= sub("_","",t2)]
  datagc[, t2:= sub(".hsq","",t2)]

  for ( i in 1:nrow(datagc)){
    f = datagc[i,file]
    tmp = fread(paste0(path,"/", f), fill=TRUE)
    datagc$cor_g[i] = tmp[Source == "rG",Variance]
    datagc$se_cor_g[i] = tmp[Source == "rG",SE]
    datagc$h2_t1[i] = tmp[Source == "V(G)/Vp_tr1",Variance]
    datagc$h2_t2[i] = tmp[Source == "V(G)/Vp_tr2",Variance]
  }
  datagc[abs(cor_g) == 1, cor_g := NA]
  datagc[abs(cor_g) == 0, cor_g := NA]

  alldatagc <- datagc[CJ(t1 = tt, t2 = tt), on = .(t1 = t1, t2 = t2)]
  # alldatagc[!is.na(cor_g), .N] / nrow(alldatagc)

  mat_corg = dcast(alldatagc, t1 ~ t2, value.var = "cor_g", fill=NA)
  setDT(mat_corg)
  setcolorder(mat_corg, tt)

  t1 = mat_corg$t1
  mat_corg[, t1:=NULL]
  mat_corg = as.matrix(mat_corg)
  rownames(mat_corg) = t1

  mat_corg = mat_corg[order(match(t1, tt)),]
  # check if mat_corg is a upper-triangle matrix
  mat_corg[lower.tri(mat_corg)] <- (t(mat_corg)[lower.tri(mat_corg)])
  # diag(mat_corg) = datah2$h2


  #### get phenotypic correlations ####
  mu <- pheno[, lapply(.SD, mean, na.rm = TRUE), by = Gruppe, .SDcols = tt]
  adjpheno <- merge(pheno, mu, by = "Gruppe", suffixes = c("", "_grpmean"))
  for (v in tt) {
    adjpheno[, (v) := get(v) - get(paste0(v, "_grpmean")) ]
  }
  # check
  x <- sample(adjpheno$kml, 1)
  y <- sample(tt, 1)
  g <- pheno[kml == x, Gruppe]
  if (! pheno[kml == x, ..y] - mu[Gruppe == g, ..y] == adjpheno[kml == x, ..y] ) stop("adj pheno seems wrong!")
  #
  tab <- adjpheno[, ..tt]
  mat_corp = cor(tab, use = "pairwise.complete.obs")
  ptab = cor.mtest(tab)
  mat_corppval = ptab$p
  # mat_corp[1:7,1:7]

  #### combine into one matrix ####
  # lower = phenotypic cor
  # upper = genetic cor
  identical(rownames(mat_corg), rownames(mat_corp))
  identical(colnames(mat_corg), colnames(mat_corp))
  if (!all(dim(mat_corg) == dim(mat_corp))) {
    stop("mat_corg and mat_corp must have the same dimensions!")
  }
  mat_corg <- as.matrix(mat_corg)
  mode(mat_corg) <- "numeric"
  mat_corp <- as.matrix(mat_corp)
  mode(mat_corp) <- "numeric"
  # define order of variables, cor and pval, BEFORE combination!
  mat_corg = mat_corg[tt, tt]
  mat_corp = mat_corp[tt, tt]
  mat_corppval = mat_corppval[tt, tt]
  # combine
  matcor = matrix(NA, nrow = nrow(mat_corg), ncol = ncol(mat_corg))
  matcor[upper.tri(matcor)] <- mat_corg[upper.tri(mat_corg)]
  matcor[lower.tri(matcor)] <- mat_corp[lower.tri(mat_corp)]
  rownames(matcor) = rownames(mat_corg)
  colnames(matcor) = colnames(mat_corg)
  diag(matcor) = datah2$h2
  # matcor[1:7,1:7]

  # all sig pval for cor_g (upper) to show
  mat_corppval[upper.tri(matcor)] = 0.01 # <0.05
  # mat_corppval[upper.tri(matcor)] <- ifelse(is.na(matcor[upper.tri(matcor)]), 1, 0.01)

  # check
  identical(rownames(matcor), rownames(mat_corppval))
  identical(colnames(matcor), colnames(mat_corppval))
  range(matcor, na.rm = TRUE)  
  sum(is.na(matcor))
  mat_corppval[which(is.na(matcor))] = NA # pval need to be NA too

  matcor = as.matrix(matcor)
  mode(matcor) <- "numeric"
  str(matcor)

  nt <- ncol(matcor)
  if (nt <5) {
    tl.cex <- 1.5  # Larger text for variable names
    cl.cex <- 1.2  # Adjust color bar labels size
    number.cex <- 1.5
    width = 4000
    height = 4000
  } else {
    tl.cex <- 1  # Default text size
    cl.cex <- 0.8  # Default color bar label size
    number.cex <- 1
    width = 4000
    height = 3800
  }
  png(args[2], width=width, height=height, res=300)  
  corrplot(matcor, p.mat = mat_corppval, sig.level=0.05, 
          method='square', addCoef.col='grey1', insig='blank', diag=TRUE, na.label = "NA",
          tl.cex = tl.cex, cl.cex = cl.cex, number.cex = number.cex, number.font = 1,
          type='full',  is.corr = FALSE)
  # method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")
  # mtext(paste0("rg_",folder,"(upper) & rp(lower)"), cex = 2, line = 1)
  dev.off()

