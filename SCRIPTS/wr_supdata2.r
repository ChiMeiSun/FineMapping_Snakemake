# args <- c("RESULTS/supdata/supdata2.txt", 
# "RESULTS/BFMAP/BW32/h2_BW32_all.txt",
# "RESULTS/BFMAP/EW40/h2_EW40_all.txt")

args <- commandArgs(TRUE)
args

library(data.table)
library(ggplot2)

files <- args[2:length(args)]
alldat <- data.table()
for (f in files){
    dat <- fread(f)
    alldat <- rbind(alldat, dat)
}

write.table(alldat, args[1], quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

png("RESULTS/supdata/supdata2.png", width = 5000, height = 3000, res = 300)
ggplot(alldat, aes(x = tragen, y = h2)) + geom_boxplot()
dev.off()