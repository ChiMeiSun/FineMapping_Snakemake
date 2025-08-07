args <- commandArgs(TRUE)
args

library(data.table)

out <- args[length(args)]
idx <- length(args) - 1
mafs <- list()

pdf(out, width=5,height=5)
for (i in 1:idx){
    gen <- sub(".frq", "", unlist(strsplit(args[i], "/"))[4])
    dat <- fread(args[i])
    boxplot(dat$MAF ~ dat$CHR, main = gen)
}
dev.off()
