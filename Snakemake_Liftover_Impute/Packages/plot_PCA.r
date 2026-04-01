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
var_explained <- eigenval/sum(eigenval)*100
head(var_explained)


id <- fread(args[3])
table(id$Generation)

data <- eigenvec[id, on = .(V1 = id)]
selcols <- c("AnimalID","id","pc1","pc2","pc3","pc4","pc5","pc6")
colnames(data)[1:8] <- selcols

data <- data[, c(selcols, "Generation"), with = FALSE][, Generation := as.factor(Generation)]

selegen <- c("Araucana_ref", "Araucana_gt", "WhiteLayer_ref", "WhiteLeghorn_ref", "Founder", "FounderBC1", "FounderBC2")
subid <- id[Generation %in% selegen, id]
subdata <- data[AnimalID %in% subid]

# col_breed = c("Araucana"="seagreen4", "Founder" = "orangered","F1"="pink2","Founder_BC1"="darkorange", "BC1"="chocolate4", 
# "Founder_BC2"="darkolivegreen3", "BC2"="azure4", "IC"="antiquewhite3")


# plot
# pdf("test.pdf", width = 7, height = 5)
pdf(args[4], width = 7, height = 5)

plot(var_explained,
     type = 'b',
     cex = 0.5,
     xlab = 'PC',
     ylab = 'variance explained [%]')

ggplot(data, aes(x = pc1, y = pc2))+
  geom_point(aes(color = Generation, shape = Generation), alpha = 0.5, size = 2)+
  geom_vline(xintercept = 0,linetype = 2, alpha = 0.7)+
  geom_hline(yintercept = 0,linetype = 2, alpha = 0.7)+
  labs(x = paste0('PC1 [',round(var_explained[1],2),'%]'),
       y = paste0('PC2 [',round(var_explained[2],2),'%]')) +
#   scale_color_manual(values = col_breed) +
  scale_shape_manual(values = 1:18) +
  theme_minimal()

ggplot(subdata, aes(x = pc1, y = pc2))+
  geom_point(aes(color = Generation, shape = Generation), alpha = 0.5, size = 2)+
  geom_vline(xintercept = 0,linetype = 2, alpha = 0.7)+
  geom_hline(yintercept = 0,linetype = 2, alpha = 0.7)+
  labs(x = paste0('PC1 [',round(var_explained[1],2),'%]'),
       y = paste0('PC2 [',round(var_explained[2],2),'%]')) +
#   scale_color_manual(values = col_breed) +
  scale_shape_manual(values = 1:18) +
  theme_minimal()

ggplot(data, aes(x = pc3, y = pc4))+
  geom_point(aes(color = Generation, shape = Generation), alpha = 0.5, size = 2)+
  geom_vline(xintercept = 0,linetype = 2, alpha = 0.7)+
  geom_hline(yintercept = 0,linetype = 2, alpha = 0.7)+
  labs(x = paste0('PC3 [',round(var_explained[3],2),'%]'),
       y = paste0('PC4 [',round(var_explained[4],2),'%]')) +
#   scale_color_manual(values = col_breed) +
  scale_shape_manual(values = 1:18) +
  theme_minimal()

ggplot(data, aes(x = pc5, y = pc6))+
  geom_point(aes(color = Generation, shape = Generation), alpha = 0.5, size = 2)+
  geom_vline(xintercept = 0,linetype = 2, alpha = 0.7)+
  geom_hline(yintercept = 0,linetype = 2, alpha = 0.7)+
  labs(x = paste0('PC5 [',round(var_explained[5],2),'%]'),
       y = paste0('PC6 [',round(var_explained[6],2),'%]')) +
#   scale_color_manual(values = col_breed) +
  scale_shape_manual(values = 1:18) +
  theme_minimal()


dev.off()



###########
# load("../DATA/geno/2024-04-19_IMAGEcomplete_CR_SCM_GGA6_Ref0_Alt1.RData")
# datgt <- data.table(PED[, c("Patient_ID", "Generation", "Breed")])
# colnames(datgt)[1] <- "id"
# datgt[Breed == "Araucana", Generation := "Araucana_gt"]
# table(datgt$Generation)
# datgt[, Breed := NULL]

# ids <- eigenvec$V1

# datref <- fread("Data/chicken_genomes_ID_Breed.csv")
# colnames(datref) <- c("id", "Generation")
# datref <- datref[id %in% ids][Generation == "Araucana", Generation := "Araucana_ref"]
# datref[grep("Layer", Generation), Generation := "WhiteLayer_ref"]
# datref[grep("Leghorn", Generation), Generation := "WhiteLeghorn_ref"]
# table(datref$Generation)

# dat <- rbind(datref, datgt)

# write.table(dat, "Data/animals_info.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')
###########
