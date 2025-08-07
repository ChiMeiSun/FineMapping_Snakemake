# {input.rdata} {output.out} {wildcards.gen0}
args <- commandArgs(TRUE)
args
load(args[1])
table(PED$Generation)
idsnp <- colnames(IMAGE)
length(idsnp) - length(PED$Patient_ID)

if (args[3] == "all") {
    write.table(cbind("IM",PED[PED$Patient_ID %in% idsnp,"Patient_ID"]),
    args[2],
    quote = FALSE, row.names = FALSE, col.names = FALSE)

} else if (args[3] == "offspring") {
    tem = PED[(PED$Generation == "F1" | PED$Generation == "BC1" | PED$Generation == "BC2" | PED$Generation == "IC"),"Patient_ID"]
    tem = tem[which(tem %in% idsnp)]
    write.table(cbind("IM",tem), args[2], quote = FALSE, row.names = FALSE, col.names = FALSE)
} else if (args[3] == "Araucana") {
    tem = PED[(PED$Breed == "Araucana"),"Patient_ID"]
    tem = tem[which(tem %in% idsnp)]
    write.table(cbind("IM",tem), args[2], quote = FALSE, row.names = FALSE, col.names = FALSE)
} else if (args[3] == "WLfounders") {
    tem = PED[(PED$Generation == "Founder"), "Patient_ID"]
    tem = tem[-c(1:6)]
    tem = tem[which(tem %in% idsnp)]
    write.table(cbind("IM",tem), args[2], quote = FALSE, row.names = FALSE, col.names = FALSE)
} else {
    write.table(cbind("IM",PED[(PED$Generation == args[3] & PED$Patient_ID %in% idsnp),"Patient_ID"]),
    args[2],
    quote = FALSE, row.names = FALSE, col.names = FALSE)

}
