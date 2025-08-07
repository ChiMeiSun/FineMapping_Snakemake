
# args = "PLOTS/manplot/combine/arrimp/BW32.png"
# args = "PLOTS/manplot/combine/arrimp_cov/BW32covEW30.png"

args <- commandArgs(TRUE)
args

# install.packages("magick")
library(magick)
library(data.table)

# PLOTS/manplot/mlma/{gen}/{pheno}_{gen}.png
# PLOTS/manplot/fastGWA/{gen}/{pheno}_{gen}.png
path = "PLOTS/manplot/"
group = unlist(strsplit(args[1],"/"))[4]
group
gens = c("all", "offspring", "IC", "BC2", "BC1", "FounderBC2", "FounderBC1")

if (group == "arrimp"){
    phe = unlist(strsplit(args[1],"/"))[5]
    phe = sub(".png", "", phe)

    files_arr = c()
    files_imp = c()
    for (g in gens){
        f1 = paste0(path,"mlma/",g,"/",phe,"_",g,".png")
        f2 = paste0(path,"fastGWA/",g,"/",phe,"_",g,".png")
        files_arr = c(files_arr, f1)
        files_imp = c(files_imp, f2)
    }
}

# PLOTS/manplot/mlma_cov/{gen}/{pheno}_{gen}.png
# PLOTS/manplot/fastGWA_cov/{gen}/{pheno}_{gen}.png
if (group == "arrimp_cov"){
    phe = unlist(strsplit(args[1],"/"))[5]
    phe = sub(".png", "", phe)

    files_arr = c()
    files_imp = c()
    for (g in gens){
        f1 = paste0(path,"mlma_cov/",g,"/",phe,"_",g,".png")
        f2 = paste0(path,"fastGWA_cov/",g,"/",phe,"_",g,".png")
        files_arr = c(files_arr, f1)
        files_imp = c(files_imp, f2)
    }
}

# read in figures
list_arr = list()
list_imp = list()

for (i in 1:length(files_arr)){
    filepath = files_arr[i]
    if (file.exists(filepath)){
            img = try(image_read(filepath))
            info = image_info(img)
            wid = info$width
            hei = info$height / 2
            img = image_crop(img, geometry = sprintf("%dx%d+0+0", wid, hei)) # crop off qqplot below
            list_arr[[i]] = img
    }
}

for (i in 1:length(files_imp)){
    filepath = files_imp[i]
    if (file.exists(filepath)){
            img = try(image_read(filepath))
            info = image_info(img)
            wid = info$width
            hei = info$height / 2
            img = image_crop(img, geometry = sprintf("%dx%d+0+0", wid, hei)) # crop off qqplot below
            list_imp[[i]] = img
    }
}

# join images
images_arr = image_join(list_arr)
images_imp = image_join(list_imp)

combined_images_arr <- image_append(images_arr, stack = FALSE)
combined_images_imp <- image_append(images_imp, stack = FALSE)

rows = ceiling( (length(files_arr)-1) / 2)
combined_images_arr <- image_montage(images_arr, tile = paste0('x',rows), geometry='3800x')
combined_images_imp <- image_montage(images_imp, tile = paste0('x',rows), geometry='3800x')
comb <- image_append(c(combined_images_arr,combined_images_imp), stack = FALSE)

montage_info <- image_info(comb)
subtitle_image <- image_blank(width = montage_info$width, height = montage_info$height/20, color = 'white')

comb <- image_append(c(subtitle_image,comb), stack = TRUE)
image_write(comb, path = args[1], format = "png")

