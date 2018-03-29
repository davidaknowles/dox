
require(data.table)

require(stringr)
require(qvalue)

require(foreach)
require(dplyr)
require(tidyr)
require(leafcutter)
require(magrittr)

source("utils.R")
source("load_data.R")

gtex_file = commandArgs(trailingOnly=TRUE)

#gtex_file="~/scailscratch/dox/GTEx/Brain_Cerebellum_100k.txt.gz"

add_loc=function(dat) {
  str_split_fixed( dat$variant_id, "_", 5 )[,1:2] %>% 
    as.data.frame(stringsAsFactors=F) %>% 
    set_colnames(c("chr","pos")) %>%
    mutate(pos=as.integer(pos)) %>% 
    cbind(dat)
}

GDRIVE_LOCATION=Sys.getenv("GDRIVE_LOCATION")

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
    rename(p_geno=`p-value`)
)

geno_threshold=1e-5

foreach(eqtl_name=names(eqtls), .combine = bind_rows) %do% {
  eqtl=eqtls[[eqtl_name]]
  eqtl_join_gtex = eqtl %>% inner_join(gtex_eqtls, by=c(gene="gene_id", RSID="RSID"))
  rep_p = eqtl_join_gtex %>% filter(p_geno < geno_threshold) %>% .$pval_nominal
  data.frame( gtex_file=gtex_file, eqtls=eqtl_name, prop_shared=nrow(eqtl_join_gtex)/nrow(eqtl), naive_rep=mean(rep_p < 0.05), p1= 1. - pi0est(rep_p)$pi0, stringsAsFactors = F)
} %>% write.table(paste0(gtex_file,"_pi0.txt"), sep="\t",row.names = F, col.names = T, quote = F)


