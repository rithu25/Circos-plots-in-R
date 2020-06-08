### Group3

library(RCircos)
setwd("~/Desktop/Circos")

####Ideogram
data("UCSC.HG38.Human.CytoBandIdeogram")
cyto.info <-UCSC.HG38.Human.CytoBandIdeogram
# print(cyto.info)
tracks.inside <-20
tracks.outside <-0
chr.exclude <-NULL
colorScales <- RCircos.Get.Heatmap.Color.Scale(heatmap.color="GreenBlackRed")

RCircos.Set.Core.Components ( cyto.info, chr.exclude, tracks.inside,tracks.outside)
#pdf("Rcircos.pdf", height = 8, width = 8, compress = T)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

### Gene layer
my_data <-read.csv("gtff3.csv")
name.col <-4
side <-"in"
track.num <-1
RCircos.Gene.Connector.Plot(genomic.data= my_data, track.num = track.num,  side=side)
track.num <-2
RCircos.Gene.Name.Plot(my_data, name.col,track.num, side)
#dev.off()


### Motifs layer
my_data2 <-read.csv("group3ciiider.csv")

my_data2$CNV <- floor(runif(n=11025, min=1, max=7))
data.col <- 6# the CNV column

name.col <- 5
track.num <- 4
side <- "in"
RCircos.Histogram.Plot(my_data2, data.col, 6, side)
#RCircos.Histogram.Plot(my_data2, data.col, track.num, side, is.sorted=FALSE, min.value=-2)
#RCircos.Line.Plot(my_data2, data.col,+ track.num, side);
#RCircos.Scatter.Plot(my_data2, data.col, track.num, side, is.sorted=FALSE)

RCircos.Gene.Connector.Plot(genomic.data= my_data2, track.num=7,  side=side)
RCircos.Gene.Name.Plot(my_data2, name.col, 8, side)
RCircos.Get.Gene.Name.Plot.Parameters()
track.num <- 6


### TAD or no TAD
data(RCircos.Scatter.Data);
my_data4 = read.csv("m3-tad-classification.csv")
data.col <- 11;
track.num <- 10;
side <- "in";
by.fold <- 0.5;
# RCircos.Scatter.Plot(my_data4, data.col,track.num, side, by.fold);

# RCircos.Line.Plot(my_data4, data.col,track.num, side);

# RCircos.Get.Plot.Colors(my_data4, color)

#### Link data

#data(RCircos.Link.Data)

RCircos.Link.Data = read.csv("lo3.csv")
RCircos.Link.Plot(RCircos.Link.Data, track.num=15, by.chromosome=TRUE)


### Log FC

my_data3 <- read.csv("log3.csv")


#$logFC <- runif(28947, -3, 3)

data.col <- 5
track.num <- 7
side <- "in"

# "by.fold" is a zero or positive number. If it's positive, then any data point with a value >= by.fold will be plotted as red color; any data point with a value <= -by.fold will be plotted as blue color; otherwise, data point will be plotted in black color.
by.fold <- 1.5

# plot scatter plot
logfc <- read.csv("log3.csv")
# logfc <- read.csv("testlog.csv")
log_loc <- logfc[,1:4]
log_ne <- logfc[,5]
log_el <- logfc[,6]
log_ne <- cbind(log_loc, log_ne)
log_el <- cbind(log_loc, log_el)

# print(log_loc[,3] - log_loc[,2])
# data("RCircos.Heatmap.Data")

# RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col=5, track.num=4, side="in")
# print(head(hmap))
# hmap[,4] <- make.unique(as.character(hmap[,4]))
# print(head(x))
# print(head(log_ne))
# print(head(log_el))
# print(log_ne)
# print(log_el)
RCircos.Heatmap.Plot(log_ne, data.col=5, track.num=15, side="in")
RCircos.Heatmap.Plot(log_el, data.col=5, track.num=16, side="in")

# tad_pos <- read.csv("tad_pos.tsv", sep="\t")[,1:3]
# tad_neg <- read.csv("tad_neg.tsv", sep="\t")[,1:3]
tad_pos_neg <- read.csv("tad_pos_neg.tsv", sep="\t")
# print(head(tad_pos_neg))
RC.param <- RCircos.Get.Plot.Parameters()
RC.param["heatmap.color"] <- "GreenBlackRed"
RCircos.Reset.Plot.Parameters(RC.param)
RCircos.Heatmap.Plot(tad_pos_neg, data.col=5, track.num=20, side="in", is.sorted=FALSE)

# RCircos.Tile.Plot(tad_pos, track.num=18, side="in", is.sorted=FALSE)
# RCircos.Tile.Plot(tad_neg, track.num=19, side="in", is.sorted=FALSE)
# RCircos.Heatmap.Plot(logfc, data.col=5, track.num=2, side="out")
# RCircos.Histogram.Plot(log_ne, data.col=5, track.num=11, side="in")
# RCircos.Histogram.Plot(log_el, data.col=5, track.num=17, side="in")
# RCircos.Line.Plot(log_ne, data.col=5, track.num=11, side="in")
# RCircos.Line.Plot(log_el, data.col=5, track.num=14, side="in")
# RCircos.Scatter.Plot(log_ne, data.col=5, track.num=11, side="in")
# RCircos.Scatter.Plot(log_el, data.col=5, track.num=17, side="in")

# data(RCircos.Line.Data, data.col=5, track.num=9, side="in")
# RCircos.Line.Plot(RCircos.Line.Data, data.col=5, track.num=9, side="in")
# RCircos.Line.Plot(log_ne, data.col=4, track.num=7, side="in")
# RCircos.Line.Plot(log_el, data.col=4, track.num=8, side="in")



#RCircos.Histogram.Plot(my_data3,data.col, track.num, side)
# RCircos.Heatmap.Plot(my_data3, data.col,
#                     track.num, side)
# RCircos.Heatmap.Plot(my_data3, data.col=5,
#                     track.num=4, side="in")

#RCircos.Parallel.Line.Plot(line.data, track.num=5, side="in")
#RCircos.Parallel.Line.Plot(line.data, line.width=2,
#                           inside.pos=2, outside.pos=2.5)
