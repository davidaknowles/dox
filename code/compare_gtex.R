
require(data.table)

require(stringr)
require(qvalue)

require(foreach)
require(dplyr)
require(tidyr)
#require(leafcutter)
require(magrittr)

source("utils.R")
source("load_data.R")

source("joint_storey_pi0.R")

gtex_file=if (interactive()) "/scratch/users/dak33/data/gtex/100k/Heart_Atrial_Appendage.allpairs.txt.gz" else commandArgs(trailingOnly=TRUE)
#gtex_file="/scratch/users/dak33/data/gtex/100k/Uterus.allpairs.txt.gz"

add_loc=function(dat) {
  str_split_fixed( dat$variant_id, "_", 5 )[,1:2] %>% 
    as.data.frame(stringsAsFactors=F) %>% 
    set_colnames(c("chr","pos")) %>%
    mutate(pos=as.integer(pos)) %>% 
    cbind(dat)
}

#GDRIVE_LOCATION=Sys.getenv("GDRIVE_LOCATION")
GDRIVE_LOCATION="/scratch/users/dak33/data/"

hg19_snps=fread(paste0("zcat < ",GDRIVE_LOCATION,"/genomes/human/snp146_maf0p05.txt.gz"), sep="\t", data.table = F) %>% 
  set_colnames(c("Ch","BP","RSID")) %>%
  filter(Ch %in% paste0("chr",1:22)) %>%
  mutate(Ch=substr(Ch,4,nchar(Ch)) %>% as.integer(), BP=BP+1)

gtex_eqtls = fread(paste0("zcat < ",gtex_file), data.table = F) %>% 
  add_loc() %>% 
  filter(chr != "X")  %>% 
  mutate(chr=as.integer(chr)) %>% 
  inner_join(hg19_snps, by=c(chr="Ch",pos="BP"))  %>% 
  mutate( gene_id=str_split_fixed(gene_id, '[.]', 2)[,1] )

#DATADIR="~/gdrive/dox_data/"
#DATADIR=paste0(gdrive_location,"dox_data/")


eqtls=list( my_eqtl = read_qtls(paste0(DATADIR,"/panama_qq_boot_1e+05/"))  %>% 
              left_join(snploc, by=c(cis_snp="snpid")) ,
  matrix_eqtl = fread(paste0("zcat < ",DATADIR,"/results0_PC10.txt.gz"), data.table = F) %>% 
    select(-FDR) %>% left_join(snploc %>% select(-pos), by=c("SNP"="snpid")) %>% 
    select(-SNP, -beta, -`t-stat`) %>% 
    rename(p_geno=`p-value`),
  ase_eqtl = fread(paste0("zcat <",DATADIR,"all_eqtl_w_ase.txt.gz"), data.table = F)  %>%
    filter(!is.na(RSID))
)

results = foreach(eqtl_name=names(eqtls), .combine = bind_rows) %do% {
  print(eqtl_name)
  eqtl=eqtls[[eqtl_name]]
  eqtl_join_gtex = eqtl %>% inner_join(gtex_eqtls, by=c(gene="gene_id", RSID="RSID")) %>%
    filter(!is.na(p_geno), !is.na(pval_nominal))
  
  print("p value thresholding")
  p_based=foreach(geno_threshold=c(1e-3, 5e-4, 1e-4, 5e-5, 1e-5), .combine = bind_rows) %do% {
    rep_p = eqtl_join_gtex %>% filter(p_geno < geno_threshold) %>% .$pval_nominal
    
    reverse_p = eqtl_join_gtex %>% filter(pval_nominal < 1e-3) %>% .$p_geno
    
    data.frame( type="p", 
                geno_threshold=geno_threshold, 
                gtex_file=gtex_file, 
                eqtls=eqtl_name, 
                prop_shared=nrow(eqtl_join_gtex)/nrow(eqtl), 
                naive_rep=mean(rep_p < 0.05), 
                naive_reverse=mean(reverse_p<0.05), 
                p1= 1. - my_pi0est(rep_p)$pi0, 
                p1_reverse=1. - my_pi0est(reverse_p)$pi0, 
                stringsAsFactors = F)
  }
  
  print("FDR thresholding")
  with_fdr=eqtl_join_gtex %>% mutate(p=p_geno, cis_snp=RSID) %>% 
    bonferroni() %>% 
    rename(RSID=cis_snp) %>%
    inner_join(gtex_eqtls, by=c(gene="gene_id", RSID="RSID"))
  q_based=foreach(geno_threshold=c(0.1, 0.05, 0.01), .combine = bind_rows) %do% {
    rep_p = with_fdr %>% filter(q < geno_threshold) %>% .$pval_nominal
    data.frame( type="q", 
                geno_threshold=geno_threshold, 
                gtex_file=gtex_file, 
                eqtls=eqtl_name, 
                prop_shared=nrow(eqtl_join_gtex)/nrow(eqtl), 
                naive_rep=NA, naive_reverse=NA, 
                p1= 1. - my_pi0est(rep_p)$pi0, 
                p1_reverse=NA, 
                stringsAsFactors = F)
  }
  
  print("Storey type approach")
  storey_joint=eqtl_join_gtex %>% filter(!is.na(p_geno), !is.na(pval_nominal)) %>% 
    select(p_geno, pval_nominal) %>% 
    joint_pi0_estimator()
  pi0=storey_joint$fit$par$pi0
  joint=data.frame( type="joint", 
              geno_threshold=NA, 
              gtex_file=gtex_file, 
              eqtls=eqtl_name, 
              prop_shared=nrow(eqtl_join_gtex)/nrow(eqtl), 
              naive_rep=pi0[2], 
              naive_reverse=pi0[3], 
              p1=storey_joint$jaccard, 
              p1_reverse=pi0[4], 
              stringsAsFactors = F)
  print("Done")
  bind_rows(p_based, q_based, joint)
}

results %>% write.table(paste0(gtex_file,"_pi0_extended.txt"), sep="\t",row.names = F, col.names = T, quote = F)


