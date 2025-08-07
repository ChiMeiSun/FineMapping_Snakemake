# args = c(mlma,manplot, pheno_gen)
# args = c("RESULTS/mlma/mlma_IC/KGW_IC.mlma","plot","KGW","IC")
args <- commandArgs(TRUE)

library(data.table)
mlma = fread(args[1], header=TRUE)
dim(mlma)
input = "Array"
if (length(grep("metaGWAS", args[1])) >0){input = "metaGWAS"}

if (length(grep("fastGWA", args[1])) >0){
     mlma = mlma[,-c("N")]
     colnames(mlma) = c("Chr", "SNP","bp","A1","A2","Freq","b","se","p")
     input = "Imputed"
}

head(mlma)

ph = args[3]
gen = args[4]

#### Manhattan plot ####
# plot blocks / SNPs
dat = data.frame(max_bp = tapply(mlma$bp,mlma$Chr, max))
cum = c(0,cumsum(dat$max_bp))
dat = cbind(Chr = rownames(dat),dat, bp_add = cum[-length(cum)])
#mlma$Chr = as.character(mlma$Chr)
dat$Chr = as.integer(dat$Chr)
mlma = merge(mlma,dat,by="Chr")
mlma$pos <- ifelse(mlma$Chr==1, mlma$bp, mlma$bp+mlma$bp_add)
dat = cbind(dat, m_bp = tapply(mlma$pos,mlma$Chr, mean)) 
# to get pos for x-axis labels (chr)

xmax = ceiling(max(mlma$pos, na.rm=TRUE) * 1.03)
xmin = floor(min(mlma$pos, na.rm=TRUE) * -0.03)


col = c('grey','black',"orangered")
col_index = rep(rep(1:2,times=100)[1:length(dat$Chr)], times = unlist(table(sort(mlma$Chr))))
# repeat 1 num.snps.chr1 times, 2 num.snps.chr2 times, 1..chr3...etc
bonf = 0.05/dim(mlma)[1]
col_index[which(mlma$p < bonf)] = 3
dat[dat$Chr == 40, "Chr"] <- "Z"    

# save plots
# png(args[2], width = 5500, height = 3500, res = 300)
# # tiff('testKGW.tiff',width = 17, height = 8, res=600, compression='lzw',units='in')
# par(mfrow = c(2, 3))

png(args[2], width = 3800, height = 5000, res = 400)
par(mfrow = c(2, 1))

par(mar = c(6, 6, 6, 2) + 0.1) # b,l,t,r

ymax = round(max(-log10(mlma$p), na.rm=TRUE)) + 2
if (is.infinite(ymax)){ymax = 50}
plot(mlma$pos, -log10(mlma$p),
     main = paste0(input,"_",ph,"_",gen), cex.main = 2,
     col = col[col_index],
     pch = 20, cex = 0.85,
     xaxt = "n", xlab = "", #rm xaxis ticks
     ylab = expression(-log[10](p)), las = 1, #horizontal y-tick-label
     cex.lab = 1.5, cex.axis = 1.5,
     xlim = c(xmin, xmax), 
     ylim = c(0, ymax),
     bty = 'n') # border type none 
axis(1, at = c(0,dat$m_bp), labels = c('', dat$Chr), 
     lty = 1, pos = -0.1, col=NA, col.ticks = 'black', cex.axis = 1.5)
abline(h = -log10(bonf), col = "steelblue") # Bonferroni correction p-val 0.05

#### qqplot ####
plot(x = -log10(ppoints(sum(!is.na(mlma$p)))),
     y = -log10(sort(mlma$p, decreasing = FALSE)),
     pch = 20,
     xlab = "Expected -log10(p)",
     ylab = "Observed -log10(p)",
     cex.lab = 1.5, cex.axis = 1.5)
abline(0,1, col = "red2")

dev.off()

