# snakemake --dag | dot -Tsvg > dag.svg
# snakemake --rulegraph | dot -Tsvg > rulegraph.svg
# snakemake --filegraph | dot -Tsvg > filegraph.svg

include: "Snakefiles/Snakefile_preps"
include: "Snakefiles/Snakefile_GWAS"
include: "Snakefiles/Snakefile_GP"
include: "Snakefiles/Snakefile_Finemap"
include: "Snakefiles/Snakefile_publish"

phenos = ['EN1','EN2','EN3','EN4','EN5','EN6','EN7','EN8','EN10','EN11','EN12','EN13','BW32','EW30','EW40','EW50','EW70']
gen0 = ['all', 'Founder', 'F1', 'FounderBC1','BC1','FounderBC2','BC2','IC','offspring','Araucana','WLfounders']
gens2 = ['all','offspring','FounderBC1','BC1','FounderBC2','BC2','IC']

wildcard_constraints:
    gen0  = "|".join(gen0)



rule all:
    input:
        "PLOTS/supdata/pheno_descriptive.pdf",
        "RESULTS/publish/pheno_descriptive.txt",
        "PLOTS/publish/PhenoCor.png",
        expand("PLOTS/publish/fig{num}.png", num=['3','4','5']),
        "PLOTS/manplot/metaGWAS/ENs_all.png",
        expand("PLOTS/manplot/mlma/{gen}/{pheno}_{gen}.png", pheno=phenos, gen = gens2),
        expand("PLOTS/manplot/fastGWA/{gen}/{pheno}_{gen}.png", pheno=phenos,gen = gens2),
        expand("PLOTS/manplot/mlma_cov/{gen}/{pheno}cov{covt}_{gen}.png", pheno=['EW30','EW40','EW50','EW70'], covt='BW32',gen=gens2),
        expand("PLOTS/manplot/fastGWA_cov/{gen}/{pheno}cov{covt}_{gen}.png", pheno=['EW30','EW40','EW50','EW70'], covt='BW32',gen=gens2),
        expand("PLOTS/manplot/mlma_cov/{gen}/{pheno}cov{covt}_{gen}.png", pheno='BW32',covt=['EW30','EW40','EW50','EW70'],gen=gens2),
        expand("PLOTS/manplot/fastGWA_cov/{gen}/{pheno}cov{covt}_{gen}.png", pheno='BW32',covt=['EW30','EW40','EW50','EW70'],gen=gens2),
        expand("PLOTS/manplot/combine/arrimp/{pheno}.png", pheno=phenos),
        expand("PLOTS/manplot/combine/arrimp_cov/{pheno}cov{covt}.png", pheno=['EW30','EW40','EW50','EW70'], covt='BW32'),
        expand("PLOTS/manplot/combine/arrimp_cov/{pheno}cov{covt}.png", pheno='BW32', covt=['EW30','EW40','EW50','EW70']),
        expand("RESULTS/gcta/sigSNP/mlma/{gen}/{pheno}_{gen}.txt", pheno=phenos, gen = gens2),
        expand("RESULTS/gcta/sigSNP/fastGWA/{gen}/{pheno}_{gen}.txt", pheno=phenos,gen = gens2),
        expand("PLOTS/BFMAP/{pheno}/{pheno}_{gen}.pdf", pheno=phenos, gen='all'),
        expand("PLOTS/BFMAP/{pheno}/recalc_{pheno}_{gen}.pdf", pheno=phenos, gen='all'),
        expand("RESULTS/supdata/SupplementaryData{num}.txt", num=['1','2']),
        expand("RESULTS/supdata/SupplementaryData{num}.xlsx", num=['3','4','5','6']),


