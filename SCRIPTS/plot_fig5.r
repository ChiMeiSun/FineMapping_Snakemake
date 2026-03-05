# args = c("PLOTS/manplot/metaGWAS/ENlate_all.png", "PLOTS/publish/fig5.png")

args = commandArgs(TRUE)
args

library(magick)
library(data.table)

# Fig. 4 metaGWAS of EN3-EN13
files = args
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

info <- image_info(images)
wid = info$width
hei = info$height
titleh = hei / 7
crop <- image_crop(images, geometry = sprintf("%dx%f+0+%f", wid, hei - titleh*1.5, titleh))

image_write(crop, path = args[length(args)], format = "png")

