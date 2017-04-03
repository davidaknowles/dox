#require(dplyr)
require(doMC)

if (F) {
    input <- read.delim("../data/gene-counts-round-two.txt.gz")
    anno <- input %>% dplyr::select(filename, individual, flow_cell, lane, index, conc)
    colnames(anno)[2]="sampleid"
    colnames(anno)[3]="individual"
    anno$filename=paste0("s",anno$filename)
    anno$sample=with(anno, paste(individual, conc, sep="_") )

    write.table(anno, file="../data/sample-anno.txt", quote=F, row.names=F, col.names = T, sep="\t")
} else {
    anno=read.table("../data/sample-anno.txt", header=T, sep="\t", stringsAsFactors = F)
}

basebase="/afs/cs.stanford.edu/u/davidknowles/scailscratch/dox/"
basedir=paste0(basebase, "bam/")

## files_exist=foreach(f=anno$filename, .combine = c) %do% {
##   file.exists(paste0(basedir,f,".bam"))
## }

foreach(sample=unique(anno$sample)) %do% {
    files=paste0(basedir,anno$filename[ anno$sample==sample ],".bam")
    outfile=paste0(basebase, "full_pileups/",sample,".pileup.gz")
    if (!file.exists(outfile)) {
        command_string=paste0("qsub -d . -q daglab -v INFILEA=",files[1],",INFILEB=",files[2],",OUTFILE=",outfile," -N ",sample," pileup.sh")
        cat(command_string,"\n")
        system( command_string )
        }
}

## fs=foreach(sample=unique(anno$sample), .combine = rbind) %do% {
##     files=paste0(basedir,anno$filename[ anno$sample==sample ],".bam")
##     outfile=paste0(basebase, "pileups/",sample,".pileup.gz")
##     c(file.size(files[1]), file.size(files[2]), file.size(outfile))
## } 

## colnames(fs)=c("bam1","bam2","pileup")
## fs=as.data.frame(fs)
## fs$bam=fs$bam1+fs$bam2
## plot(fs$bam, fs$pileup)
