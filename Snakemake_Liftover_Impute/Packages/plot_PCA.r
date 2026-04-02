# args <- c(
#      "Output/PCA/merged_pruned.eigenvec",
#      "Output/PCA/merged_pruned.eigenval",
#      "Data/animals_info.txt",
#      "Output/PCA/PCA.pdf"
# )

args = commandArgs(TRUE)
args

library(data.table)
library(ggplot2)

eigenvec <- fread(args[1])
eigenval <- unlist(fread(args[2]))
eigenval[eigenval < 0] <- 0
var_explained <- eigenval/sum(eigenval)*100
head(var_explained)


id <- fread(args[3])
table(id$Generation)

data <- eigenvec[id, on = .(V1 = id)]
selcols <- c("AnimalID","id","pc1","pc2","pc3","pc4","pc5","pc6")
colnames(data)[1:8] <- selcols

data <- data[, c(selcols, "Generation"), with = FALSE]

# selegen <- c("Araucana_ref", "Araucana_gt", "WhiteLayer_ref", "WhiteLeghorn_ref", "Founder", "FounderBC1", "FounderBC2")
# subid <- id[Generation %in% selegen, id]
# subdata <- data[AnimalID %in% subid]

ref <- c("Araucana_ref", "WhiteLayer_ref", "WhiteLeghorn_ref")
target <- setdiff(unique(data$Generation), ref)
data$Generation <- factor(data$Generation, levels = c(ref, target))
table(data$Generation)

color_palette <- c(
  # Reference populations 
  "Araucana_ref" = "gray50",     
  "WhiteLeghorn_ref" = "gray60",  
  "WhiteLayer_ref" = "gray70",    
  
  # Target populations 
  "Araucana" = "seagreen4",      
  "Founder" = "orangered",  
  "F1" = "pink2",         
  "FounderBC1" = "darkorange",  
  "BC1" = "chocolate4",        
  "FounderBC2" = "darkolivegreen3",  
  "BC2" = "violet",        
  "IC" = "antiquewhite3"          
)

# plot
# pdf("test.pdf", width = 7, height = 5)
pdf(args[4], width = 7, height = 5)

plot(var_explained,
     type = 'b',
     cex = 0.5,
     xlab = 'Principal Component',
     ylab = 'variance explained (%)')

ggplot() +
  geom_point(data = data, 
    aes(x = pc1, y = pc2, color = Generation, shape = Generation), 
    alpha = 0.5, size = 2) +
  scale_color_manual(values = color_palette, name = "Reference/Generation") +
  scale_shape_manual(values = 1:20, guide = "none") +
  guides(
    color = guide_legend(
      override.aes = list(shape = 1:length(unique(data$Generation))))  
  ) +
  geom_vline(xintercept = 0,linetype = 2, alpha = 0.7)+
  geom_hline(yintercept = 0,linetype = 2, alpha = 0.7)+
  labs(x = paste0('PC1 [',round(var_explained[1],2),'%]'),
       y = paste0('PC2 [',round(var_explained[2],2),'%]')) +
  theme_minimal()

# ggplot(subdata, aes(x = pc1, y = pc2))+
#   geom_point(aes(color = Generation, shape = Generation), alpha = 0.5, size = 2)+
#   geom_vline(xintercept = 0,linetype = 2, alpha = 0.7)+
#   geom_hline(yintercept = 0,linetype = 2, alpha = 0.7)+
#   labs(x = paste0('PC1 [',round(var_explained[1],2),'%]'),
#        y = paste0('PC2 [',round(var_explained[2],2),'%]')) +
# #   scale_color_manual(values = col_breed) +
#   scale_shape_manual(values = 1:18) +
#   theme_minimal()

ggplot() +
  geom_point(data = data, 
    aes(x = pc3, y = pc4, color = Generation, shape = Generation), 
    alpha = 0.5, size = 2) +
  scale_color_manual(values = color_palette, name = "Reference/Generation") +
  scale_shape_manual(values = 1:20, guide = "none") +
  guides(
    color = guide_legend(
      override.aes = list(shape = 1:length(unique(data$Generation))))  
  ) +
  geom_vline(xintercept = 0,linetype = 2, alpha = 0.7) +
  geom_hline(yintercept = 0,linetype = 2, alpha = 0.7) +
  labs(x = paste0('PC3 [',round(var_explained[3],2),'%]'),
       y = paste0('PC4 [',round(var_explained[4],2),'%]')) +
#   scale_color_manual(values = col_breed) +
  theme_minimal()

# ggplot(data, aes(x = pc5, y = pc6))+
#   geom_point(aes(color = Generation, shape = Generation), alpha = 0.5, size = 2)+
#   geom_vline(xintercept = 0,linetype = 2, alpha = 0.7)+
#   geom_hline(yintercept = 0,linetype = 2, alpha = 0.7)+
#   labs(x = paste0('PC5 [',round(var_explained[5],2),'%]'),
#        y = paste0('PC6 [',round(var_explained[6],2),'%]')) +
# #   scale_color_manual(values = col_breed) +
#   scale_shape_manual(values = 1:18) +
#   theme_minimal()

dev.off()



###########
# load("../DATA/geno/2024-04-19_IMAGEcomplete_CR_SCM_GGA6_Ref0_Alt1.RData")
# datgt <- data.table(PED[, c("Patient_ID", "Generation", "Breed")])
# colnames(datgt)[1] <- "id"
# datgt[Breed == "Araucana", Generation := "Araucana"]
# table(datgt$Generation)
# datgt[, Breed := NULL]


# # SAMEA6873070	IM003
# # SAMEA6873071	IM005
# # SAMEA6873072	IM006
# # SAMEA6873073	IM007
# # SAMEA6873074	IM008
# # SAMEA6873075	IM009
# ids <- eigenvec$V1
# ids_dupARU <- c("SAMEA6873070","SAMEA6873071","SAMEA6873072","SAMEA6873073","SAMEA6873074","SAMEA6873075")
# ids <- setdiff(ids, ids_dupARU)

# datref <- fread("Data/chicken_genomes_ID_Breed.csv")

# colnames(datref) <- c("id", "Generation")
# datref <- datref[id %in% ids][Generation == "Araucana", Generation := "Araucana_ref"]
# datref[grep("Layer", Generation), Generation := "WhiteLayer_ref"]
# datref[grep("Leghorn", Generation), Generation := "WhiteLeghorn_ref"]
# table(datref$Generation)

# dat <- rbind(datref, datgt)

# write.table(dat, "Data/animals_info.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')
# write.table(dat[, .(id,id)], "Data/pca_samples.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
###########
