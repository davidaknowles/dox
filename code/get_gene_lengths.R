
library(GenomicRanges)
library(rtracklayer)
#library(Rsamtools)

#Load the annotation and reduce it
GTF <- import.gff("~/Dropbox/splicing/leafcutter/leafcutter/data/gencode.v26.annotation.gtf.gz",
                  format="gtf", genome="GRCm38.71", feature.type="exon") #, asRangedData=F)
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))

elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
get_width <- function(x) { sum(elementMetadata(x)$widths) }
output <- sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), get_width)

write.table(output, file="gene_lengths.tsv", sep="\t", quote=F, col.names = F)
