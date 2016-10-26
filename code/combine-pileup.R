# Pull allelic counts from the .pileup.gz files generated via submit-mpileup.R and combine into
# a single file. 

require(doMC)
registerDoMC(7)
require(data.table)

anno=read.table("../data/sample-anno.txt", header=T, sep="\t", stringsAsFactors = F)

basebase="/afs/cs.stanford.edu/u/davidknowles/scailscratch/dox/"
basedir=paste0(basebase, "bam/")

samples=unique(anno$sample)

all_hets=foreach(sample=samples, .combine = c, .errorhandling = "remove") %dopar% {
    outfile=paste0(basebase, "pileups/",sample,".pileup.gz")
    if (!file.exists(outfile)) return(NULL)
    b=fread(paste0("zcat < ",outfile))
    setDF(b)
    colnames(b)=c("chr","pos","r","y","ref","alt")
    paste0(b$chr,":",b$pos)
}

ta=table(all_hets)
to_keep=names(ta)[ta>1]

all_samp=rbindlist( foreach(sample=samples, .errorhandling = "remove") %dopar% {
    outfile=paste0(basebase, "pileups/",sample,".pileup.gz")
    if (!file.exists(outfile)) return(NULL)
    b=fread(paste0("zcat < ",outfile))
    setDF(b)
    colnames(b)=c("chr","pos","r","y","ref","alt")
    b$sample=sample
    b[ paste0(b$chr,":",b$pos) %in% to_keep, ]
} )
setDF(all_samp)

gzf=gzfile("../data/ase.txt.gz", "w")
write.table(all_samp, file=gzf, quote=F, row.names=F, col.names = T, sep="\t")
close(gzf)
