# args = c("RESULTS/PCA/impute/bcftools_lifted_imputed_gt_QC_pruned.eigenvec", "RESULTS/PCA/impute/bcftools_lifted_imputed_gt_QC_pruned.eigenval",
# "DATA/geno/2024-04-19_IMAGEcomplete_CR_SCM_GGA6_Ref0_Alt1.RData",
# "PLOTS/PCA/bcf_imputed_IMAGE_PCA.pdf")


args <- commandArgs(TRUE)
args


library(data.table)
library(ggplot2)

eigenvec <- fread(args[1])
eigenval <- unlist(fread(args[2]))
var_explained <- eigenval/sum(eigenval)*100
var_explained


load(args[3])
PED = as.data.table(PED)
colnames(PED)
PED$Generation[1:6] = "Araucana"
eigenvec = eigenvec[PED[,c(1,14)], on = .(V2 = Patient_ID)]
colnames(eigenvec)[1:8] = c("famid","sid","pc1","pc2","pc3","pc4","pc5","pc6")
colnames(eigenvec)[ncol(eigenvec)] = "breed"
col_breed = c("Araucana"="seagreen4", "Founder" = "orangered","F1"="pink2","Founder_BC1"="darkorange", "BC1"="chocolate4", 
"Founder_BC2"="darkolivegreen3", "BC2"="azure4", "IC"="antiquewhite3")


pdf(args[4], height=5, width=7)

plot(var_explained,
     type = 'b',
     cex = 0.5,
     xlab = 'PC',
     ylab = 'variance explained [%]')

ggplot(eigenvec)+
  geom_point(aes(x = pc1, y = pc2, col = breed))+
  geom_vline(xintercept = 0,linetype = 2)+
  geom_hline(yintercept = 0,linetype = 2)+
  labs(x = paste0('PC1 [',round(var_explained[1],2),'%]'),
       y = paste0('PC2 [',round(var_explained[2],2),'%]')) +
  scale_color_manual(values = col_breed)

ggplot(eigenvec)+
  geom_point(aes(x = pc3, y = pc4, col = breed))+
  geom_vline(xintercept = 0,linetype = 2)+
  geom_hline(yintercept = 0,linetype = 2)+
  labs(x = paste0('PC3 [',round(var_explained[3],2),'%]'),
       y = paste0('PC4 [',round(var_explained[4],2),'%]')) +
  scale_color_manual(values = col_breed)

ggplot(eigenvec)+
  geom_point(aes(x = pc5, y = pc6, col = breed))+
  geom_vline(xintercept = 0,linetype = 2)+
  geom_hline(yintercept = 0,linetype = 2)+
  labs(x = paste0('PC5 [',round(var_explained[5],2),'%]'),
       y = paste0('PC6 [',round(var_explained[6],2),'%]')) +
  scale_color_manual(values = col_breed)

dev.off()
