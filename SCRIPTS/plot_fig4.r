# args = c("PLOTS/manplot/mlma/all/BW32_all.png", "PLOTS/manplot/mlma/all/EW40_all.png",
# "PLOTS/manplot/fastGWA/all/BW32_all.png", "PLOTS/manplot/fastGWA/all/EW40_all.png", 
# "PLOTS/manplot/fastGWA_cov/all/BW32covEW40_all.png","PLOTS/manplot/fastGWA_cov/all/EW40covBW32_all.png",
# "PLOTS/publish/fig4.png")

args = commandArgs(TRUE)
args

library(magick)
library(data.table)

# Fig. 4 Array and impute GWAS for BW32 and EW40, group “all” + impute GWAS with covariate for BW32 and EW40, group “all”
files = args[c(1,3,5,2,4,6,7)]

files = files[-length(files)]
rows = ceiling( (length(files)) / 3)
list = list()
for (i in 1:length(files)){
    filepath = files[i]
    if (file.exists(filepath)){
        img = try(image_read(filepath))
        info = image_info(img)
        wid = info$width
        hei = info$height / 2
        img = image_crop(img, geometry = sprintf("%dx%d+0+0", wid, hei)) # crop off qqplot below
        list[[i]] = img
    }
} 

images = image_join(list)
combined_image <- image_append(images, stack = FALSE)
combined_image <- image_montage(images, tile = paste0('x',rows), geometry='3800x') # tile=col x row

montage_info <- image_info(combined_image)
montage_info
subtitle_image <- image_blank(width = montage_info$width, height = montage_info$height/20, color = 'white')

comb <- image_append(c(subtitle_image,combined_image), stack = TRUE)
image_write(comb, path = args[length(args)], format = "png")

