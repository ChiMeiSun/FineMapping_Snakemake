# args = c("RESULTS/gcta/mlma/all/BW32_all.mlma","PLOTS/manplot/mlma/all/BW32_all.png","BW32","all")

args <- commandArgs(TRUE)
args

library(data.table)

mlmap <- args[1]
outplot <- args[2]
ph <- args[3]
gen <- args[4]

mlma <- fread(mlmap, header = TRUE)
dim(mlma)

main_txt <- paste0("Array_",ph,"_",gen)
if (length(grep("imp", mlmap)) > 0) main_txt <- paste0("Imputed_",ph,"_",gen)
if (length(grep("metaGWAS", mlmap)) > 0) main_txt <- paste0("metaGWAS_EN1-EN13_",gen)


#### Manhattan plot ####
chr_max <- mlma[, .(max_bp = max(bp)), by = Chr]
setorder(chr_max, Chr)
chr_max[, bp_add := c(0, cumsum(max_bp[-.N]))]

mlma <- merge(mlma, chr_max, by = "Chr")
mlma[, pos := bp + bp_add]

chr_means <- mlma[, .(m_bp = mean(pos)), by = Chr]

# get plot limits
xmax <- ceiling(max(mlma$pos, na.rm = TRUE) * 1.03)
xmin <- floor(min(mlma$pos, na.rm = TRUE) * -0.03)
ymax <- round(max(-log10(mlma$p), na.rm=TRUE)) + 2
if (is.infinite(ymax)) ymax <- 50

bonf <- 0.05 / nrow(mlma)

mlma[, col_index := as.integer(factor(Chr)) %% 2 + 1]  # 1 for odd, 2 for even
mlma[p < bonf, col_index := 3]
col_index <- mlma$col_index
col = c('grey','black',"orangered")

chr_means[, Chr := as.character(Chr)]
chr_means[Chr == "40", Chr := "Z"]


# save plots
png(outplot, width = 9.5, height = 12, res = 300, units = "in")
# tiff(outplot, width = 9, height = 13, res = 300, units = "in") # much larger file than png
par(mfrow = c(2, 1))
par(mar = c(6, 6, 6, 2) + 0.1) # b,l,t,r

plot(mlma$pos, -log10(mlma$p),
     main = main_txt, cex.main = 2,
     col = col[col_index],
     pch = 20, cex = 0.85,
     xaxt = "n", xlab = "", # rm xaxis ticks
     ylab = expression(-log[10](p)), las = 1, #horizontal y-tick-label
     cex.lab = 1.5, cex.axis = 1.5,
     xlim = c(xmin, xmax), 
     ylim = c(0, ymax),
     bty = 'n') # border type none 
axis(1, at = c(0, chr_means$m_bp), labels = c('', chr_means$Chr), 
     lty = 1, pos = -0.1, col = NA, col.ticks = 'black', cex.axis = 1.5)
abline(h = -log10(bonf), col = "steelblue", lwd = 2) # Bonferroni correction p-val 0.05

#### qqplot ####
p_valid <- mlma[!is.na(p), p]

plot(x = -log10(ppoints(length(p_valid))),
     y = -log10(sort(p_valid)),
     pch = 20,
     xlab = "Expected -log10(p)",
     ylab = "Observed -log10(p)",
     cex.lab = 1.5, cex.axis = 1.5)
abline(0, 1, col = "red2")

dev.off()

