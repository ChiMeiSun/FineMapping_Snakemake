cd IMAGE_Publish/SCRIPTS/packages/bfmap_0.65/example

../bfmap --compute_grm 1 --binary_genotype_file geno --snp_info_file all.snp_info.csv --output grm1 --num_threads 10

../bfmap --varcomp --phenotype phen.csv --trait milk --binary_grm_file grm1 --output milk --num_threads 10

(../bfmap --sss --phenotype phen.csv --trait milk --snp_info_file topQTL.snp_info.csv --snp_weight weight --binary_genotype_file geno --binary_grm grm1 --heritability 0.307879 --output milk.topQTL --num_threads 10) &> sss.log

# GWAS
../bfmap --assoc --phenotype phen.csv --trait milk --snp_info_file all.snp_info.csv --snp_set single --binary_genotype_file geno --binary_grm grm1 --heritability 0.307879 --output milk.single.assoc --num_threads 10

# test annot_file, cat_prob file
R --vanilla
library(data.table)
snp = fread("milk.topQTL.pip.csv")
cons = c("missense_variant","intron_variant","intergenic_variant")
annot = data.table(SNPname = snp$SNPname, 
                    cons = cons[sample(1:3, 200, replace = TRUE)])

catp = data.table(cons = cons, prior = 1, freq = 1)

write.table(annot, "test_annot.txt", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(catp, "test_catprob.txt", quote=FALSE, col.names=TRUE, row.names=FALSE)
q()