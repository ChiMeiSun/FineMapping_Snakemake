# args = c("RESULTS/reml/IMAGE/BW32_EW40.hsq","RESULTS/reml/files/commands.txt")
args = commandArgs(TRUE)
args

library(data.table)
path = sub(".hsq","", args[1])
cmds = fread(args[2], sep="\n", header=FALSE)

ind = grep(path, cmds$V1)
cmd = as.character(cmds[ind,])
system(cmd)


