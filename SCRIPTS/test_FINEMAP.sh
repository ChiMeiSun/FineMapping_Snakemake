# wget -P SCRIPTS/packages http://www.christianbenner.com/finemap_v1.4.2_x86_64.tgz
# wget -P SCRIPTS/packages http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz
# tar -xzvf SCRIPTS/packages/ldstore_v2.0_x86_64.tgz -C SCRIPTS/packages/
# SCRIPTS/packages/finemap_v1.4.2_x86_64/finemap
# 4:72624837-77912432
        # prepare GWAS statistics of interested region
            awk 'BEGIN {{FS=OFS=" "}} $1==4 && $3>=72624837 && $3<=77912432 {{print $2,$1,$3,$4,$5,$7,$8,$9,$10}}' \
            RESULTS/gcta/mlma_impute/all/EW50_all.mlma | sort -k9,9n > test/1.z
            awk '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9=0}}' test/1.z > test/2.z
            echo 'MAF>0.5'
            awk '$6>0.5 {{print}}' test/2.z | wc -l

            echo 'rsid chromosome position allele1 allele2 maf beta se flip' > test/_modi.z            
            awk '{{if ($6>0.5) {{$6=1-$6; t=$4; $4=$5; $5=t; $9=1}} print}}' test/2.z >> test/_modi.z  
            echo num_snp_flip
            awk '$9==1 {{print}}' test/_modi.z | grep -v rsid | wc -l

        # calculate SNPs pearson's correlation (using same SNPs and samples as GWAS)
            awk 'NR>1 {{print $1}}' test/_modi.z > test/_snplist.txt

            awk '{{print $2}}' DATA/pheno/EW50.txt | sort > test/EW50.txt
            awk '{{print $2}}' INTERMEDIATE/IDlist/IDlist_all.txt | sort > all.txt
            comm -12 test/EW50.txt all.txt > test/EW50_all.txt
            awk '{{print "IM", $0}}' test/EW50_all.txt > test/EW50_all2.txt
            echo num_samples
            wc -l test/EW50_all2.txt

            plink --bfile INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC --r --matrix --out test/ \
            --extract test/_snplist.txt \
            --keep test/EW50_all2.txt --chr-set 40 --allow-extra-chr
            
            echo num_snps_from_mlma
            wc -l test/_snplist.txt
            echo num_snps_plink_extract_LD
            wc -l test/.ld

            ns=512
            echo 'z;ld;snp;config;cred;log;n_samples' > test/_sss
            echo 'test/_modi.z;test/.ld;test/_sss.snp;test/_sss.config;test/_sss.cred;test/_sss.log;'$ns \
            >> test/_sss

            SCRIPTS/packages/finemap_v1.4.2_x86_64/finemap --sss --n-threads 5 --log \
            --in-files BW32_all

            echo 'z;ld;snp;config;cred;log;n_samples' > test/_cond
            echo 'test/_modi.z;test/.ld;test/_cond.snp;test/_cond.config;test/_cond.cred;test/_cond.log;'$ns \
            >> test/_cond
            SCRIPTS/packages/finemap_v1.4.2_x86_64/finemap --cond --flip-beta \
            --in-files test/_cond

            leadsnp=$(awk 'NR==2 {{print $2}}' test/_sss.snp)
            echo $leadsnp
            line=$(awk -v pattern=$leadsnp '$0 == pattern {{print NR}}' test/_snplist.txt)
            echo $line
            awk -v line_number=$line 'NR==line_number {{print $0}}' test/.ld > test/_leadsnp.ld
            
            rm test/.ld
            rm test/1.z
            rm test/2.z
            rm test/_modi.z


        awk '{{$1=$2; print}}' INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC.fam > tmp1.txt
        mv tmp1.txt INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC.fam


plink2 --bfile "INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC" --export bgen-1.2 --out "INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC"  --chr-set 40 --allow-extra-chr --threads 2
bgenix -g "INTERMEDIATE/beagleimpute/bcftools_lifted_imputed_gt_QC.bgen" -index

SCRIPTS/packages/ldstore_v2.0_x86_64/ldstore --in-files RESULTS/FINEMAP/BW32/all/masterLD --write-bcor --read-only-bgen --n-threads 3
SCRIPTS/packages/ldstore_v2.0_x86_64/ldstore --in-files RESULTS/FINEMAP/BW32/all/masterLD --bcor-to-text --n-threads 3
SCRIPTS/packages/finemap_v1.4.2_x86_64/finemap --sss --n-threads 5 --log --in-files RESULTS/FINEMAP/BW32/all/masterFM
SCRIPTS/packages/finemap_v1.4.2_x86_64/finemap --cond --log --in-files RESULTS/FINEMAP/BW32/all/masterFM