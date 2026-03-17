# wget -P SCRIPTS/packages http://www.christianbenner.com/finemap_v1.4.2_x86_64.tgz
# wget -P SCRIPTS/packages http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz
# tar -xzvf SCRIPTS/packages/ldstore_v2.0_x86_64.tgz -C SCRIPTS/packages/
# SCRIPTS/packages/finemap_v1.4.2_x86_64/finemap
# 4:72624837-77912432

plink2 --bfile "INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC" --export bgen-1.2 --out "INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC"  --chr-set 40 --allow-extra-chr --threads 2
bgenix -g "INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC.bgen" -index

SCRIPTS/packages/ldstore_v2.0_x86_64/ldstore --in-files RESULTS/FINEMAP/BW32/all/masterLD --write-bcor --read-only-bgen --n-threads 3
SCRIPTS/packages/ldstore_v2.0_x86_64/ldstore --in-files RESULTS/FINEMAP/BW32/all/masterLD --bcor-to-text --n-threads 3
SCRIPTS/packages/finemap_v1.4.2_x86_64/finemap --sss --n-threads 5 --log --in-files RESULTS/FINEMAP/BW32/all/masterFM
SCRIPTS/packages/finemap_v1.4.2_x86_64/finemap --cond --log --in-files RESULTS/FINEMAP/BW32/all/masterFM


# BFMAP
# wget -P SCRIPTS/packages https://github.com/jiang18/bfmap/releases/download/v0.91/bfmap-v0.91-x86_64-linux.zip
# unzip SCRIPTS/packages/bfmap-v0.91-x86_64-linux.zip -d SCRIPTS/packages 

## create snpinfo file
awk -v chrom=4 -v start=74175424 -v end=77608972 '
    $1 == chrom && $4 >= start && $4 <= end {print $2}
' INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC.bim > tmp_4:74175424-77608972.txt

echo id > tmp_4:74175424-77608972.csv
cat tmp_4:74175424-77608972.txt >> tmp_4:74175424-77608972.csv
head tmp_4:74175424-77608972.csv

## write tmp_pheno.csv
Rscript --vanilla SCRIPTS/prep_bfmap.r DATA/pheno/Phenotypes_Final.csv DATA/geno/2024-04-19_IMAGEcomplete_CR_SCM_GGA6_Ref0_Alt1.RData


## make grm (QTL genotypes for h2)
SCRIPTS/packages/bfmap_0.91/bfmap --compute_grm 2 --binary_genotype_file INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC             --snp_info_file tmp_4:74175424-77608972.csv             --output_file tmp_4:74175424-77608972


## run BFMAP

## regional QTL h2 estimation
SCRIPTS/packages/bfmap_0.91/bfmap --varcomp --phenotype_file tmp_pheno.csv --trait_name EW70             --binary_grm_file tmp_4:74175424-77608972             --covariate_file INTERMEDIATE/covar/covar_bfmap.csv             --output RESULTS/BFMAPsss/EW70/VC/4:74175424-77608972_EW70_all

H2=$(grep "proportion" RESULTS/BFMAPsss/EW70/VC/4:74175424-77608972_EW70_all.varcomp.csv | cut -d',' -f2)
echo heritability: $H2


# forward selection for fine mapping
SCRIPTS/packages/bfmap_0.91/bfmap --sss --phenotype_file tmp_pheno.csv --trait_name EW70             --binary_genotype_file INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC             --snp_info_file tmp_4:74175424-77608972.csv             --binary_grm_file INTERMEDIATE/BFMAPsss/GRM/allsnp_all             --heritability $H2             --covariate_file INTERMEDIATE/covar/covar_bfmap.csv             --output RESULTS/BFMAPsss/EW70/FS/4:74175424-77608972_EW70_all