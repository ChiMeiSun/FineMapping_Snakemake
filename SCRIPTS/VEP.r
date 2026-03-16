# args = c("RESULTS/BFMAP/EW40/credsets_EW40_all.txt", 
# "Snakemake_Liftover_Impute/Output/Imputed/bcftools_lifted_imputed_gt_DR2.vcf.gz",
# "RESULTS/BFMAP/EW40/vep_EW40_all.vcf", "EW40", "all")

args = commandArgs(TRUE)
args

library(data.table)
#
credp <- args[1]
vcfp <- args[2]
vepp <- args[3]
pheno <- args[4]
gen <- args[5]
#

cred <- fread(credp)
tragen = paste0(pheno,"_",gen)

system(paste0("rm -r temp_",tragen))
system(paste0("mkdir temp_",tragen))

region <- cred[, .(Chr, Pos)]
regpath <- paste0("temp_",tragen,"/region.txt")
write.table(region, regpath, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

regvcfpath <- paste0("temp_",tragen,"/region.vcf")
system(
    paste0("bcftools view -R ",regpath," ",vcfp," -Ov -o ",regvcfpath)
    ,wait=TRUE
)
print("Finished bcftools subset region.vcf")

param = "--species gallus_gallus_gca000002315v5 --cache --dir_cache /data/users/chi.sun/.vep --merged --force_overwrite --everything --pick --stats_text --stats_html"
system(
    paste0("vep -i ",regvcfpath," -o ",vepp," ",param)
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