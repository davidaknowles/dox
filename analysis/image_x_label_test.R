
library(ggplot2)  ## devtools::install_github("hadley/ggplot2)
library(grid)     ## rasterGrob
library(EBImage)  ## readImage
library(ggthemes) ## theme_minimal

## ##########
## INDEPENDENT CODE TO BE SOURCED:
## ##########
# user-level interface to the element grob
my_axis = function(img) {
  structure(
    list(img=img),
    class = c("element_custom","element_blank", "element") # inheritance test workaround
  )
}
# returns a gTree with two children: the text label, and a rasterGrob below
element_grob.element_custom <- function(element, x,...)  {
  stopifnot(length(x) == length(element$img))
  tag <- names(element$img)
  # add vertical padding to leave space
  g1 <- textGrob(paste0(tag, "\n\n\n\n\n"), x=x, vjust=0.6)
  g2 <- mapply(rasterGrob, x=x, image=element$img[tag], 
               MoreArgs=list(vjust=0.7, interpolate=FALSE,
                             height=unit(3,"lines")),
               SIMPLIFY=FALSE)
  
  gTree(children=do.call(gList, c(g2, list(g1))), cl="custom_axis")
}
# gTrees don't know their size and ggplot would squash it, so give it room
grobHeight.custom_axis = heightDetails.custom_axis = function(x, ...)
  unit(6, "lines")
## ##########
## END
## ##########

## ##########
## OBTAIN FLAGS:
## ##########
library(rvest)

## GDP per capita, top 10 countries
url      <- "https://en.wikipedia.org/wiki/List_of_countries_by_GDP_(nominal)_per_capita"
html     <- read_html(url)
gdppc    <- html_table(html_nodes(html, "table")[3])[[1]][1:10,]

## clean up; remove non-ASCII and perform type conversions
gdppc$Country <- gsub("Â ", "", gdppc$Country)
gdppc$Rank    <- iconv(gdppc$Rank,    "latin1", "ASCII", sub="")
gdppc$Country <- iconv(gdppc$Country, "latin1", "ASCII", sub="")
gdppc$`US$`   <- as.integer(sub(",", "", gdppc$`US$`))

## flag images (yes, this processing could be done neater, I'm sure)
## get the 200px versions
flags_img  <- html_nodes(html_nodes(html, "table")[3][[1]], "img")[1:10]
flags_url  <- paste0('http://', sub('[0-9]*px', '200px', sub('\\".*$', '', sub('^.*src=\\"//', '', flags_img))))
flags_name <- sub('.*(Flag_of)', '\\1', flags_url)

if(!dir.exists("flags")) dir.create("flags")
for(flag in seq_along(flags_url)) {
  switch(Sys.info()[['sysname']],
         Windows= {download.file(flags_url[flag], destfile=file.path("flags", paste0(flag,"_", flags_name[flag])), method="auto", mode="wb")},
         Linux  = {download.file(flags_url[flag], destfile=file.path("flags", paste0(flag,"_", flags_name[flag])))},
         Darwin = {download.file(flags_url[flag], destfile=file.path("flags", paste0(flag,"_", flags_name[flag])))})
}
## ##########
## END
## ##########

## load the images from filenames
npoints <- length(flags_name)

pics  <- vector(mode="list", length=npoints)
image.file <- dir("flags", full.names=TRUE)
image.file <- image.file[order(as.integer(sub("_.*","",sub("flags/","",image.file))))]
for(i in 1:npoints) {
  pics[[i]] <- EBImage::readImage(image.file[i])
}
names(pics) <- sub(".svg.png","",sub(".*Flag_of_","",image.file))

## create a dummy dataset

y       <- gdppc$`US$`
x       <- names(pics)
dat     <- data.frame(x=factor(x, levels=names(pics)), y=y)

## create the graph, as per normal now with @baptiste's adapted grob processing
## NB: #85bb65 is the color of money in the USA apparently.
gg <- ggplot(dat, aes(x=x, y=y/1e3L, group=1)) 
gg <- gg + geom_bar(col="black", fill="#85bb65", stat="identity")
gg <- gg + scale_x_discrete()
gg <- gg + theme_minimal()
gg <- gg + scale_fill_discrete(guide=FALSE)
gg <- gg + theme(plot.background = element_rect(fill="grey90"))
gg <- gg + labs(title="GDP per capita", 
                subtitle="Top 10 countries", 
                x="", y="$US/1000", 
                caption=paste0("Source: ",url))
gg <- gg + theme(axis.text.x  = my_axis(pics), ## that's much better
                 axis.text.y  = element_text(size=14),
                 axis.title.x = element_blank())
gg

