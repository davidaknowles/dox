
library("dplyr")
library("tidyr")
library(data.table)
source("utils.R")

require(rstan)
source("load_data.R")

registerDoMC(7)

require(suez)

geneloc=read.table(paste0(DATADIR,"genelocGRCh38.txt.gz"),header=T,stringsAsFactors = F)

sample_kernel=readRDS("../data/Kern.rds")

anno_for_suez= anno %>% select(-individual) %>% rename(individual=dbgap, condition=conc) %>% mutate(condition=factor(condition))

rarg_mutants= c("4714_4827", "7666_460f")

input=quantile_normalize(input)
 
genes=intersect(rownames(input),geneloc$geneid)
rownames(geneloc)=geneloc$geneid

eigen_sample_kernel=eigen(sample_kernel)

no_geno = model.matrix( ~ condition, data=anno_for_suez) # [,2:5]
  
N=ncol(input)

geno=as.numeric(anno_for_suez$individual %in% rarg_mutants)
  
results=foreach(gene=genes, .combine = bind_rows) %dopar% {
    
    y=input[gene,] %>% as.numeric
    y=y-mean(y)
    
    data=list(N=N,U_transpose_x=t(eigen_sample_kernel$vectors) %*% no_geno,P=ncol(no_geno), U_transpose_y=t(eigen_sample_kernel$vectors) %*% y %>% as.numeric, lambda=eigen_sample_kernel$values)
    
    init=list(sigma2=0.1, sigma2_k=1.0, beta=lm(y ~ no_geno - 1) %>% coef )
    
    fit_no_geno=suez:::stan_optimizing_wrapper(suez:::stanmodels$suez_step_2, data, init=init, as_vector=F)
    
      lrt = function(data) {
        data$U_transpose_x=t(eigen_sample_kernel$vectors) %*% cbind( no_geno, geno )
        data$P=ncol(data$U_transpose_x)
        init=fit_no_geno$par
        init$beta=c(init$beta,0.0)
        
        fit_geno=suez:::stan_optimizing_wrapper(suez:::stanmodels$suez_step_2, data, init=init, as_vector=F )
        
        interact=model.matrix(~geno:condition,data=anno_for_suez)
        interact=interact[,3:ncol(interact)]
        data_interact=data
        data_interact$U_transpose_x=t(eigen_sample_kernel$vectors) %*% cbind( no_geno, geno, interact )
        data_interact$P=ncol(data_interact$U_transpose_x)
        
        init=fit_geno$par
        init$beta=c(init$beta,numeric(ncol(interact)))
        fit_interact=suez:::stan_optimizing_wrapper(suez:::stanmodels$suez_step_2, data_interact, init=init, as_vector=F)
        
        list( fit_geno=fit_geno, fit_interact=fit_interact )
      }
      
      lrt_true=lrt( data )
      fit_geno=lrt_true$fit_geno
      fit_interact=lrt_true$fit_interact
      
      data.frame(gene=gene, l0=fit_no_geno$value, l_geno=fit_geno$value, l_interact=fit_interact$value, stringsAsFactors=F  )
}

results = results %>% mutate( p_geno=lrt_pvalue(l_geno-l0,df=1),
                p_interact=lrt_pvalue(l_interact-l_geno,df=4) ) %>% 
  select(-starts_with("l"))

results %<>% mutate(q_geno=p.adjust(p_geno), q_interact=p.adjust(p_interact))

results_file= gzfile( "../supp_data/rarg_trans_eqtl.txt.gz","w")
results %>% arrange( pmin(p_interact,p_geno)) %>% format(digits=5) %>% write.table(results_file, quote = F, row.names = F, col.names = T, sep="\t")
close(results_file)



 

sum(results$q_geno < 0.05)
sum(results$q_interact < 0.05)

ensg_to_hugo=read.table("../data/ensg_to_hugo.txt.gz",header=T,stringsAsFactors = F,sep="\t") %>%
  rename(hugo=Approved.Symbol, gene=Ensembl.Gene.ID)
hits=results %>% filter(q_geno < 0.05 | q_interact < 0.05) %>% left_join( ensg_to_hugo, by="gene" ) %>% arrange(q_interact)

require(gridExtra)
pdf("../figures/rarg_hits.pdf", width=10, height=8)
do.call( grid.arrange, c( foreach(i=1:nrow(hits)) %do% {
  hit=hits[i,]
  gene_name=hit$hugo
  genotype_labels=c("GG","GA","AA")
  ylabel=bquote( .(gene_name) ~ expression ~ ("log"[2]~cpm) )
  data.frame(y= input[hit$gene,] %>% as.numeric(), geno=factor(geno, 0:2,genotype_labels), conc=anno$conc) %>% filter(!is.na(geno)) %>% ggplot(aes(as.factor(conc), y, col=geno, alpha=geno))  + geom_point(position = position_jitterdodge(dodge.width =  .75, jitter.width = 0.25, jitter.height = 0.), size=3) + ylab(ylabel) + xlab( if(i!=6) expression("Doxorubicin concentration ("*mu*"M)")) + theme_bw(base_size=13) + scale_color_manual(values=cbPalette, name="rs2229774", guide=if(i==6) "legend" else F) +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom", legend.background = element_blank() ) + scale_alpha_discrete(range=c(0.5,1), guide=F)+ ggtitle(ifelse(hit$q_interact<hit$q_geno,paste0("Response eQTL p=",format(hit$q_interact,digits=1)),paste0("Marginal eQTL p=",format(hit$q_geno,digits=1))))
}, nrow=2) ) 
dev.off()

data.frame(y= input[hit$gene,] %>% as.numeric(), geno=factor(geno, 0:2,genotype_labels), conc=anno$conc) %>% filter(!is.na(geno)) %>% ggplot(aes(as.factor(conc), y, col=geno, alpha=geno))  + geom_point(position = position_jitterdodge(dodge.width =  .75, jitter.width = 0.25, jitter.height = 0.), size=3) + ylab(ylabel) + xlab(expression("Doxorubicin concentration ("*mu*"M)")) + theme_bw(base_size=13) + scale_color_manual(values=cbPalette, name="rs2229774") +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "right", legend.background = element_blank() ) + scale_alpha_discrete(range=c(0.5,1), guide=F)+ ggtitle(ifelse(hit$q_interact<hit$q_geno,paste0("Response eQTL p=",format(hit$q_interact,digits=1)),paste0("Marginal eQTL p=",format(hit$q_geno,digits=1))))
ggsave("../figures/rarg_hits_legend.pdf", width=10/3, height = 8/2)

# Don't see any relationship with conserved rarg 

# PAQR3 https://www.ncbi.nlm.nih.gov/pubmed/27212020
# https://www.nature.com/articles/s41598-017-05276-2
# LRRC2 Mitochondrial regulator? http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0170458
# SGIP1 is neuronal?! May play a role in the regulation of energy homeostasis.
# VMA21: Myopathy with excessive autophagy protein, probably not very relevant? 
# NMRK1: enzyme involved in metabolism. Bit complex? https://www.nature.com/articles/ncomms13103

