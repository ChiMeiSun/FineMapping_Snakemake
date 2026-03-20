

# args = c("Output/Phased/conformgt/IMAGE_phased_aligned.vcf.gz",
#     "Output/Phased/bcftools_lifted_qc_phased.bref3",
#     "Packages/beagle.22Jul22.46e.jar","2","5","500000","15",
#     "Output/check_imputation/CV_gt_imp_rep2_fold5.pdf",
#     "Output/check_imputation/CV_gt_imp_rep2_fold5.txt",
#     "Output/check_imputation/CV_gt_imp_rep2_fold5_MAFbin.txt")
args <- commandArgs(TRUE)
args

library(data.table)
set.seed(2024)

swap = function(x){
    res = ifelse(x==0, 2, ifelse(x==2,0, ifelse(x==1,1,1)))
    return(res)
}

gt = args[1]
vcf = fread(gt)
colnames(vcf)[1] = "chr"
snps = vcf[,c("chr","POS")]

ref = args[2]

jar = args[3]
reps = as.numeric(args[4])
folds = as.numeric(args[5])
mem = as.numeric(args[6])
mem = format(mem, scientific = FALSE, trim = TRUE)
threads = as.numeric(args[7])

para = "parallel --jobs 0 --use-cpus-instead-of-cores"
prep = paste(c(1:reps), collapse=" ")
pfold = paste(c(1:folds), collapse=" ")
param = paste0("::: ",prep," ::: ",pfold)
system("rm -r tmp")
system("mkdir tmp")

for (i in 1:reps){
    nsnpf <- floor(nrow(snps)/folds) # fsnps+- per fold
    nsnpl <- nrow(snps) %% nsnpf
    ind <- c(rep(1:folds,each = nsnpf), sample(1:folds,nsnpl) )
    ind <- sample(ind)
    print(head(ind))
    print(table(ind))

    for(k in 1:folds){
        snps_fold <- snps[which(ind==k),]
        res[rep==i & fold==k, fsnp := nrow(snps_fold)]
        table(snps_fold$chr)
        
        write.table(snps_fold, paste0("tmp/snps_rep",i,"fold",k,".txt"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
    }
    print("finish writing masked SNPs for this repetition")
}

# cmd="bcftools-1.20/bcftools view -T ^tmp/snps_rep1fold1.txt Output/Phased/conformgt/IMAGE_phased_aligned.vcf.gz -Oz -o tmp/maskfoldgt_rep1fold1.vcf.gz"
# system(cmd)
system(
    paste0(para," bcftools-1.20/bcftools view -T ^tmp/snps_rep{1}fold{2}.txt ",gt," -Oz -o tmp/maskfoldgt_rep{1}fold{2}.vcf.gz ",param)
    ,wait=TRUE)

# cmd="java -Xmx$((500000/ 1000 - 2 ))g -jar Packages/beagle.22Jul22.46e.jar gt=tmp/maskfoldgt_rep1fold1.vcf.gz ref=Output/Phased/bcftools_lifted_qc_phased.bref3 out=tmp/maskfoldgt_imp_rep1fold1"
# system(cmd)
system(
    paste0(para," java -Xmx$((",mem,"/ 1000 - 2 ))g -jar ",jar," gt=tmp/maskfoldgt_rep{1}fold{2}.vcf.gz ref=",ref," out=tmp/maskfoldgt_imp_rep{1}fold{2} ",param)
    ,wait=TRUE)

# cmd="bcftools-1.20/bcftools index -f tmp/maskfoldgt_imp_rep1fold1.vcf.gz"
# system(cmd)
system(
    paste0(para," bcftools-1.20/bcftools index -f tmp/maskfoldgt_imp_rep{1}fold{2}.vcf.gz ",param)
    ,wait=TRUE)

# cmd="bcftools-1.20/bcftools view -R tmp/snps_rep1fold1.txt tmp/maskfoldgt_imp_rep1fold1.vcf.gz -Oz -o tmp/maskfoldgt_imp_fsnp_rep1fold1.vcf.gz"
# system(cmd)
system(
    paste0(para," bcftools-1.20/bcftools view -R tmp/snps_rep{1}fold{2}.txt tmp/maskfoldgt_imp_rep{1}fold{2}.vcf.gz -Oz -o tmp/maskfoldgt_imp_fsnp_rep{1}fold{2}.vcf.gz ",param)
    ,wait=TRUE)

pdf(args[8], width=14, height=9)
par(mfrow = c(2,3))

res <- data.table(rep = rep(1:reps,each=folds),
    fold = rep(1:folds,times=reps))
res_bin <- data.table()


for (i in 1:reps){
    for(k in 1:folds){
        title = paste0("rep: ",i,", fold: ",k)
        print(title)
       
        impf <- fread(paste0("tmp/maskfoldgt_imp_fsnp_rep",i,"fold",k,".vcf.gz"))
        colnames(impf)[1]= "chr"
        if (sum(duplicated(impf$POS)) != 0) warning("Duplicated position detected")
        
        masksnps <- fread(sprintf("tmp/snps_rep%sfold%s.txt", i, k))
        orif = merge(vcf, masksnps, by.x = c("chr", "POS"), by.y = c("V1", "V2"))

        merge = orif[,c(1:9)][impf[,1:9],on=.(POS=POS, chr=chr), nomatch=NULL]

        # fsnps in imp 
        res[rep==i & fold==k, fsnps_imp := nrow(impf)]

        # imp snps not found in fsnps
        res[rep==i & fold==k, fsnps_idk := sum(!impf$ID %in% merge$i.ID)]
       
        tab = unique(merge, by="ID")    
        # snps for analysis: overlap of masked snps from ori and imp(no dups)
        res[rep==i & fold==k, fsnps_ana := nrow(tab)]
        # fsnps missing
        res[rep==i & fold==k, fsnps_miss := nrow(masksnps) - fsnps_ana]

        tab[, refimp := substr(i.REF,1,1)]
        tab[, altimp := substr(i.ALT,1,1)]
        p_match = round(nrow(tab[REF == refimp & ALT == altimp,]) / nrow(tab),3) # fsnps REF/ALT of gt/imp match
        res[rep==i & fold==k, prop_snpmatch := p_match]
        p_match
        tab[!(REF == refimp & ALT == altimp),]

        ori012 = orif[ID %in% tab$ID,-c(1:9)]
        imp012 = impf[ID %in% tab$i.ID,-c(1:9)]

        if (dim(ori012)[1] != dim(imp012)[1]) warning("Number of SNPs in fold do not match") # fsnps match?
        identical(colnames(ori012), colnames(imp012)) # samples match betw ori & imp
        ori012 = ori012[, lapply(.SD, function(column) {
          as.numeric(substr(column, 1, 1)) + as.numeric(substr(column, 3, 3))
        })]
        imp012 = imp012[, lapply(.SD, function(column) {
          as.numeric(substr(column, 1, 1)) + as.numeric(substr(column, 3, 3))
        })]

        resv = tab[,c(1:5,10:12,15,17,18)]
        for (v in 1:nrow(resv)) {
            oriv = factor(ori012[v,], level=0:2)
            impv = factor(imp012[v,], level=0:2)
            cor = round(cor(as.numeric(oriv), as.numeric(impv)),3)
            tabv = table(impv, oriv)
            p0 = sum(diag(tabv))/sum(tabv)
            tot = data.table(timp = apply(tabv,1,sum), tori = apply(tabv,2,sum))
            pc = tot[,sum(timp*tori)]/sum(tabv)^2
            iqs = (p0-pc)/(1-pc)

            DR2 = sub("DR2=","",unlist(strsplit(resv[v, i.INFO], ";"))[1])
            af = sub("AF=","",unlist(strsplit(resv[v, i.INFO], ";"))[2])
            resv[v, DR2 := ..DR2]
            resv[v, af := ..af]
            resv[v, cor := ..cor]
            resv[v, p0 := round(..p0,3)]
            resv[v, pc := round(..pc,3)]
            resv[v, iqs := round(..iqs,3)]
        }

        resv[, DR2 := as.numeric(DR2)]
        resv[, af := as.numeric(af)]
        resv[, maf := pmin(af, 1-af)]
        resv[, cut:= cut(resv$maf, breaks=c(-0.0001,0.01,0.05,0.1,0.3,0.5))]
        plot(resv$cut, resv$DR2, main=title, xlab="MAF", ylab="DR2")
        plot(resv$cut, resv$p0, main=title, xlab="MAF", ylab="Concordance")
        plot(resv$cut, resv$iqs, main=title, xlab="MAF, all masked SNPs", ylab="IQS")
        plot(resv[DR2>=0.8,cut], resv[DR2>=0.8,iqs], main=title, xlab="MAF, DR2>=0.8", ylab="IQS")
        plot(resv$cut, resv$cor, main=title, xlab="MAF", ylab="Correlation")
        plot(resv$DR2, resv$cor, main=title, xlab="DR2", ylab="Correlation")

        resv[cor<=-0.8,]
        res[rep==i & fold==k, iqs:= round(mean(resv$iqs, na.rm=TRUE),3)]
        res[rep==i & fold==k, p0:= round(mean(resv$p0),3)]
        res[rep==i & fold==k, m_cor:= round(mean(resv$cor, na.rm=TRUE),3)]
        res[rep==i & fold==k, min_cor:= round(min(resv$cor, na.rm=TRUE),3)]
        res[rep==i & fold==k, DR2:= round(mean(resv$DR2),3)]

        res_bin <- rbind(res_bin,
            resv[, .(
                iqs = round(mean(iqs, na.rm=TRUE),3),
                p0 = round(mean(p0),3),
                m_cor = round(mean(cor, na.rm=TRUE),3),
                min_cor = round(min(cor, na.rm=TRUE),3),
                DR2 = round(mean(DR2),3)
                )
                , by = cut][, `:=`(
                    rep = i,
                    fold = k)
                ][order(cut)]
        )

        print(res[rep==i & fold==k,])
    }
}
dev.off()

write.table(res, args[9], col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(res_bin, args[10], col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

