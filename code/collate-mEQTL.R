require(dplyr)
require(doMC)

anno <- read.delim("../data/sample_annotation.txt", stringsAsFactors = F)

NPCs=c(0,1,2,5,10,15,20)

concs=sort(unique(anno$conc))

basedir="~/scailscratch/dox/mEQTL_results/"

foreach (conc=concs, .combine=cbind) %do% {
    foreach(npc=NPCs,.combine = c) %do% {
        suffix=paste0(conc,"_PC",npc,".txt.gz")
        print(suffix)
        filename=paste0(basedir,"results",suffix)
        if (!file.exists(filename)) return(NA)
        qtls=read.table(filename ,header=T)
        myq=qtls %>% group_by(gene) %>% summarize( SNP=SNP[which.min(p.value)], p=min(p.value)*length(p.value))
        myq$p[myq$p>1]=1
        myq$q=p.adjust(myq$p,method="BH")
        gzf=gzfile(paste0(basedir,"summary",suffix),"w")
        write.table(as.data.frame(myq), gzf, quote=F, sep="\t", row.names=F)
        close(gzf)
        nsig=sum(myq$q<0.05,na.rm=T)
        print(nsig)
        nsig
    }
}

conc=0.0
npc=10
suffix=paste0("_pc_",conc,"_PC",npc,".txt.gz")
print(suffix)
filename=paste0(basedir,"results",suffix)
require(data.table)
qtls=fread(paste("zcat <",filename),header=T)
colnames(qtls)[5]="pvalue"
myq=qtls %>% group_by(gene) %>% summarize( SNP=SNP[which.min(pvalue)], p=min(pvalue)*length(pvalue))
myq$p[myq$p>1]=1
myq$q=p.adjust(myq$p,method="BH")
gzf=gzfile(paste0(basedir,"summary",suffix),"w")
write.table(as.data.frame(myq), gzf, quote=F, sep="\t", row.names=F)
close(gzf)
nsig=sum(myq$q<0.05,na.rm=T)
print(nsig)
nsig
