# https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM
# R script to read the GRM binary file
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}

grm_arr <- ReadGRMBin("RESULTS/gcta/GRM/GRM_bin_IMAGE_QC_all")
str(grm_arr)

grm_imp <- ReadGRMBin("RESULTS/gcta/GRM/GRM_bin_imputed_QC_all")
str(grm_imp)

print(paste0("cor_diag: ", cor(grm_arr$diag, grm_imp$diag)))
print(paste0("cor_lower_triangle: ", cor(grm_arr$off, grm_imp$off)))

# args = c("RESULTS/BFMAP/EW70/credsets_EW70_all.txt",      
# "RESULTS/BFMAP/EW70/vep_EW70_all.vcf",               
# "DATA/geno/ensembl_genes_GRCg6a.txt",                
# "RESULTS/BFMAPrecalc/EW70/credsets_EW70_all.txt",    
# "RESULTS/BFMAPrecalc/EW70/geneprobs_EW70_all.txt",   
# "RESULTS/BFMAPrecalc/EW70/enrichcat_EW70_all.txt",   
# "RESULTS/BFMAPrecalc/EW70/geneprobsori_EW70_all.txt",
# "PLOTS/BFMAP/EW70/recalc_EW70_all.pdf")

# args = c("RESULTS/supdata/SupplementaryData6.xlsx",        
# "RESULTS/BFMAPrecalc/BW32/geneprobs_BW32_all.txt",
# "RESULTS/BFMAPrecalc/EW40/geneprobs_EW40_all.txt",
# "RESULTS/BFMAPrecalc/EW50/geneprobs_EW50_all.txt",
# "RESULTS/BFMAPrecalc/EW70/geneprobs_EW70_all.txt",
# "RESULTS/BFMAPrecalc/EN2/geneprobs_EN2_all.txt",  
# "RESULTS/BFMAPrecalc/EN3/geneprobs_EN3_all.txt",  
# "RESULTS/BFMAPrecalc/EN4/geneprobs_EN4_all.txt",  
# "RESULTS/BFMAPrecalc/EN5/geneprobs_EN5_all.txt",  
# "RESULTS/BFMAPrecalc/EN6/geneprobs_EN6_all.txt",  
# "RESULTS/BFMAPrecalc/EN7/geneprobs_EN7_all.txt",  
# "RESULTS/BFMAPrecalc/EN8/geneprobs_EN8_all.txt",  
# "RESULTS/BFMAPrecalc/EN10/geneprobs_EN10_all.txt",
# "RESULTS/BFMAPrecalc/EN11/geneprobs_EN11_all.txt",
# "RESULTS/BFMAPrecalc/EN12/geneprobs_EN12_all.txt",
# "RESULTS/BFMAPrecalc/EN13/geneprobs_EN13_all.txt")

# args = c("RESULTS/supdata/SupplementaryData6.xlsx",        
# "RESULTS/BFMAPrecalc/BW32/geneprobsori_BW32_all.txt",
# "RESULTS/BFMAPrecalc/EW40/geneprobsori_EW40_all.txt",
# "RESULTS/BFMAPrecalc/EW50/geneprobsori_EW50_all.txt",
# "RESULTS/BFMAPrecalc/EW70/geneprobsori_EW70_all.txt",
# "RESULTS/BFMAPrecalc/EN2/geneprobsori_EN2_all.txt",  
# "RESULTS/BFMAPrecalc/EN3/geneprobsori_EN3_all.txt",  
# "RESULTS/BFMAPrecalc/EN4/geneprobsori_EN4_all.txt",  
# "RESULTS/BFMAPrecalc/EN5/geneprobsori_EN5_all.txt",  
# "RESULTS/BFMAPrecalc/EN6/geneprobsori_EN6_all.txt",  
# "RESULTS/BFMAPrecalc/EN7/geneprobsori_EN7_all.txt",  
# "RESULTS/BFMAPrecalc/EN8/geneprobsori_EN8_all.txt",  
# "RESULTS/BFMAPrecalc/EN10/geneprobsori_EN10_all.txt",
# "RESULTS/BFMAPrecalc/EN11/geneprobsori_EN11_all.txt",
# "RESULTS/BFMAPrecalc/EN12/geneprobsori_EN12_all.txt",
# "RESULTS/BFMAPrecalc/EN13/geneprobsori_EN13_all.txt")