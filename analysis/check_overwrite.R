
bams=read.table('~/scailscratch/dox/bam/ls.temp',sep="",header=F,skip=2)
head(bams)

genecounts=read.table('~/scailscratch/dox/counts/ls.temp',sep="",header=F)
head(genecounts)

summary(bams)
summary(genecounts)

