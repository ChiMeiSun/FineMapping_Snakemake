# args <- c("DATA/geno/2024-04-19_IMAGEcomplete_CR_SCM_GGA6_Ref0_Alt1.RData", "DATA/pheno/Phenotypes_Final.csv",
# "RESULTS/overview/pheno_descriptive.txt", "PLOTS/overview/pheno_descriptive.pdf", "RESULTS/supdata/SupplementaryData1.txt")
args <- commandArgs(TRUE)
args

library(data.table)
library(ggplot2)
library(patchwork) 

load(args[1])
pheno <- fread(args[2])

color_gen = c("FounderBC1"="darkorange", "BC1"="chocolate4", "FounderBC2"="darkolivegreen3", "BC2"="azure4", "IC"="antiquewhite3")
gens = c("FounderBC1","FounderBC2","BC1","BC2","IC")

#### pheno statstics table ####
colnames(pheno)
tab = pheno[,c(10:17,19:27)]
stat = tab[, lapply(.SD, function(x){
              N = sum(!is.na(x))
              Mean = round(mean(x, na.rm=TRUE), 2)
              SD = round(sd(x, na.rm=TRUE), 2)
              CV = round((SD / Mean) *100, 2)

              Min = min(x, na.rm=TRUE)
              Max = max(x, na.rm=TRUE)
              pct_min = mean(x == Min, na.rm=TRUE) * 100
              pct_max = mean(x == Max, na.rm=TRUE) * 100
              Min_fmt = sprintf("%.f(%.1f%%)", Min, pct_min)
              Max_fmt = sprintf("%.f(%.1f%%)", Max, pct_max)
              list(N=N,Mean=Mean,SD=SD,CV=CV,Min=Min_fmt,Max=Max_fmt)
            }), .SDcols = colnames(tab)]

statab <- rbindlist(lapply(names(stat), function(trait) {
  data.table(Trait = trait, t(unlist(stat[[trait]])))
}), fill = TRUE)
# Set column names
colnames(statab) = c("Trait","N", "Mean", "SD", "CV(%)", "Min(%)", "Max(%)")
statab$Weeks = c(paste0(seq(20,48,4),"-",seq(23,52,4)), paste0(seq(56,68,4),"-",seq(59,72,4)), 32,30,40,50,70)
setcolorder(statab, c("Trait","N", "Weeks"))

write.table(statab, args[3], quote=FALSE, col.names=TRUE, row.names=FALSE)

####   ####


colnames(PED)
colnames(PED)[5] = "kml"
pheno[, kml:=as.character(kml)]
tab = merge(pheno, PED[,c(5,14)], by = "kml")
colnames(tab)
tab = tab[,c(10:17,19:27,29)]

write.table(tab, args[5], quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
## percentage NA
perc_na <- function(x){
  round(sum(is.na(x))/length(x)*100)
}
tab_percna = tab[, lapply(.SD, perc_na), by = Generation]
tab_percna = melt(tab_percna, id.vars = "Generation", variable.name = "pheno", value.name = "perc_NA")
get_na_label = function(x) {
  return(data.frame(y = x+2, 
                    label = paste(x)))
}


pdf(args[4], width = 20, height = 8)
# Percentage of NAs for each phenotype
g1 = ggplot(tab_percna, aes(x = pheno, y = perc_NA, fill = Generation)) + 
  geom_point(shape = 21, size = 4, position = position_dodge(width = 0.5))+ 
  ggtitle("Percentage of NAs")+
  scale_fill_manual(values = color_gen)+
  stat_summary(fun.data = get_na_label, geom = "text", hjust = 0, size = 3.5) +
  theme(text = element_text(size = 16), axis.text = element_text(size = 10), 
  axis.title.x = element_text(vjust = -1)) +
  theme(plot.margin = unit(c(0.5,0.5,1,1), "cm")) # trbl
  

# BW32 by each generation
bw = tab[,c("BW32","Generation")]
bw$Generation = factor(bw$Generation, levels = gens)
ggbw= ggplot(bw, aes(x = Generation, y = BW32, fill = Generation)) +
      geom_violin(scale="width") + 
      scale_fill_manual(values = color_gen) +
      labs(y = "BW32", text = element_text(size = 16)) +
      theme(plot.margin = unit(c(0.5,0.5,1,1), "cm")) # trbl

# EWs by each generation
ew = tab[,c("EW30","EW40","EW50","EW70","Generation")]
ew$Generation = factor(ew$Generation, levels = gens)
ew = melt(ew, id.vars = "Generation",
  measure.vars = c("EW30","EW40","EW50","EW70"), variable.name="EWs")

ggew = ggplot(ew, aes(x = Generation, y = value, fill=Generation)) +
      geom_violin(scale="width") + 
      scale_fill_manual(values = color_gen) +
      facet_grid(EWs~.) +
      labs(y = "EW(g)", text = element_text(size = 16)) +
      theme(plot.margin = unit(c(0.5,0.5,1,1), "cm")) # trbl

# ENs by each generation
en = tab[,c(paste0("EN",c(1:8,10:13)),"Generation")]
en$Generation = factor(en$Generation, levels = gens)
en = melt(en, id.vars = "Generation",
  measure.vars = paste0("EN",c(1:8,10:13)), variable.name="ENs")

ggen = ggplot(en, aes(x = ENs, y = value, fill=Generation)) +
      geom_violin(scale="width") + 
      scale_fill_manual(values = color_gen) +
      facet_grid(Generation~.) +
      labs(y = "EN", text = element_text(size = 16)) +
      theme(plot.margin = unit(c(0.5,0.5,1,1), "cm")) # trbl

print(g1)
print(ggbw)
print(ggew)
print(ggen)

dev.off()
