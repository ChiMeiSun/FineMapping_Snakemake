# args = c("DATA/geno/2024-04-19_IMAGEcomplete_CR_SCM_GGA6_Ref0_Alt1.RData","DATA/pheno/Phenotypes_Final.csv",
# "RESULTS/GRM/GRM_bin_IMAGE_QC_all","RESULTS/GRM/GRM_bin_imputed_QC_all",
# "RESULTS/gcta/bireml/files/pheno.txt","RESULTS/gcta/bireml/files/filehsq.txt","RESULTS/gcta/bireml/files/commands.txt")
args = commandArgs(TRUE)
args
# need binary grm
# gcta --bfile INTERMEDIATE/ori/IMAGE_QC --make-grm-bin --out grm --autosome-num 40

library(data.table)
tt = c('EN1','EN2','EN3','EN4','EN5','EN6','EN7','EN8','EN10','EN11','EN12','EN13','BW32','EW30','EW40','EW50','EW70')

load(args[1])
PED = setDT(PED)
PED[,Patient_KM := as.integer(Patient_KM)]

pheno = fread(args[2])
pheno[,kml := as.integer(kml)]
pheno = pheno[PED[,c(1,5)], on=.(kml = Patient_KM), nomatch=NULL]
selecols = c("Patient_ID", tt)
pheno = pheno[, ..selecols]
pheno = cbind("IM", pheno)
write.table(pheno, args[5], row.names=FALSE, col.names=FALSE, quote=FALSE)

grmar = args[3]
grmimp = args[4]

identical(tt, colnames(pheno)[-c(1:2)])
cmd = data.table()
files = data.table()

nt = length(tt)
for ( i in 1:(nt-1) ){
    t1 = tt[i]
    for ( j in ((i+1):nt) ){
        t2 = tt[j]

        tn = paste0(t1,"_",t2)
        pathar = paste0("RESULTS/gcta/bireml/IMAGE/",tn)
        pathimp = paste0("RESULTS/gcta/bireml/Imputed/",tn)
        cmd = rbind(cmd, paste0("gcta --reml-bivar ",i," ",j," --grm ",grmar," --pheno ",args[5]," --out ",pathar))
        cmd = rbind(cmd, paste0("gcta --reml-bivar ",i," ",j," --grm ",grmimp," --pheno ",args[5]," --out ",pathimp))
        files = rbind(files, paste0(pathar,".hsq"))
        files = rbind(files, paste0(pathimp,".hsq"))
       
    }
}

write.table(files, args[6], row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(cmd, args[7], row.names=FALSE, col.names=FALSE, quote=FALSE)

