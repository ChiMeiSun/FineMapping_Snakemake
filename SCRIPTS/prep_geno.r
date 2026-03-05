# args=c({input.rdata} {input.pheno} {output.vcf})
args <- commandArgs(TRUE)
args

library(data.table)
load(args[1])
pheno <- fread(args[2])

IMAGE[1:5,1:5]
sum(is.na(IMAGE))
PED[1:10,]
IMAP6[1:5,]
identical(rownames(IMAGE),IMAP6$RS)

# make vcf
# 0 = 0/0 
# 1 = 0/1 
# 2 = 1/1 
# NA = ./.

#vcf[1:5,1:5]
geno = IMAGE
sum(is.na(geno))
lookup2 <- c("0/0","0/1", "1/1")
vv = lookup2[geno+1] 
vv[is.na(vv)] <- "./."
dim(vv) = dim(geno)
# pp2 = ifelse(ped == 0, "0 0", ifelse(ped == 1, "1 0", ifelse(ped == 2, "2 2", ifelse(ped == -1, "-1 -1", ped)))) 
# 2 minute

rownames(vv) = rownames(geno)
colnames(vv) = colnames(geno)
vv[1:5,1:5]

identical(IMAP6$RS, rownames(vv))
table(IMAP6$RA)
vcf = data.frame(
    CHROM = IMAP6$CHR, 
    POS = IMAP6$POS,
    ID = IMAP6$RS, 
    REF = substr(IMAP6$RA,1,1),
    ALT = substr(IMAP6$RA,3,3),
    QUAL = ".", FILTER = ".", INFO = ".", FORMAT = "GT")
vcf[vcf$CHROM == "Z", "CHROM"] = 40

vcf = cbind(vcf,vv)
colnames(vcf)[1] = "#CHROM"
vcf[1:5,1:15]

header = c("##fileformat=VCF",
paste0("##filedate=",Sys.Date()),
"##source=SCRIPTS/prep_geno.r",
paste0("##reference=",args[1]),
"##phasing=none",
"##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>")
writeLines(header, args[3])
write.table(vcf, args[3], append=TRUE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# # make ped: 
# # ref.ref = 0 = 1,1
# # ref.alt = 1 = 1,2
# # alt.alt = 2 = 2,2
# # NA = 0,0
# lookup <- c("1\t1","1\t2", "2\t2")
# #pp = lookup[as.character(ped)] # 1 minute
# pp = lookup[ped+1] 
# pp[is.na(pp)] <- "0\t0"
# dim(pp) = dim(ped)
# # pp2 = ifelse(ped == 0, "0 0", ifelse(ped == 1, "1 0", ifelse(ped == 2, "2 2", ifelse(ped == -1, "-1 -1", ped)))) 
# # 2 minute

# p = data.frame(
#     fam = substr(rownames(ped)[1], 1, 2),
#     id = rownames(ped), 
#     pid = PED$Father_ID[match(rownames(ped), PED$Patient_ID)], 
#     mid = PED$Mother_ID[match(rownames(ped), PED$Patient_ID)],
#     sex = PED$Sex[match(rownames(ped), PED$Patient_ID)], 
#     phe = -9)

# ped = cbind(p,pp)
# dim(ped)
# ped[1:20,1:10]

# # make map
# map = data.frame(
#     chr = IMAP6$CHR, 
#     rs = IMAP6$RS, 
#     gd = 0, 
#     pos = IMAP6$POS)
# table(IMAP6$RA)
# map[map$chr == "Z", "chr"] = 40
