# args = c("RESULTS/BFMAP/EW40/credsets_EW40_all.txt", 
# "Snakemake_Liftover_Impute/Output/Imputed/bcftools_lifted_imputed_gt_DR2.vcf.gz",
# "RESULTS/BFMAP/EW40/vep_EW40_all.vcf")

args = commandArgs(TRUE)
args

x = unlist(strsplit(args[1], "[/_]"))
x
pheno = x[3]
gen = sub(".txt", "", x[length(x)])
tragen = paste0(pheno,"_",gen)

system(paste0("rm -r temp_",tragen))
system(paste0("mkdir temp_",tragen))

library(data.table)
cred = fread(args[1])
cred
vcf = args[2]

region <- cred[, .(Chr, Pos)]
regpath <- paste0("temp_",tragen,"/region.txt")
write.table(region, regpath, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

regvcfpath <- paste0("temp_",tragen,"/region.vcf")
system(
    paste0("bcftools view -R ",regpath," ",vcf," -Ov -o ",regvcfpath)
    ,wait=TRUE
)
print("Finished bcftools subset region.vcf")

veppath <- args[3]
param = "--species gallus_gallus_gca000002315v5 --cache --dir_cache /data/users/chi.sun/.vep --merged --force_overwrite --everything --pick --stats_text --stats_html"
system(
    paste0("vep -i ",regvcfpath," -o ",veppath," ",param)
    ,wait=TRUE
)
print("vep of region.vcf finished")

system(paste0("rm -r temp_",tragen))

# vep_install --AUTO cf --SPECIES gallus_gallus_gca000002315v5_merged --DESTDIR /data/users/chi.sun/.vep 
# --DESTDIR not working
# mkdir /data/users/chi.sun/.vep/
# cp -r /home/chi.sun/.vep/* /data/users/chi.sun/.vep/
# du -sh /home/chi.sun/.vep/*
# du -sh /data/users/chi.sun/.vep/*

# 160 : gallus_gallus_gca000002315v5_merged_vep_114_GRCg6a.tar.gz (3 GB)
# 161 : gallus_gallus_gca000002315v5_refseq_vep_114_GRCg6a.tar.gz (3 GB)
# 162 : gallus_gallus_gca000002315v5_vep_114_GRCg6a.tar.gz (717 MB)