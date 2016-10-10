require(dplyr)

basedir="~/scailscratch/dox/mEQTL_results/"

args=if (interactive()) c("0.0","1") else commandArgs(trailingOnly = T)

conc=as.numeric(args[1])
npc=as.integer(args[2])

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
