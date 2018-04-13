
require(data.table)
require(foreach)
require(dplyr)
require(tidyr)
require(leafcutter)
require(magrittr)
source("../code/utils.R")
source("../code/load_data.R")

gdrive_location=Sys.getenv("GDRIVE_LOCATION")
#gdrive_location="~/gdrive/"
serie=read.csv(paste0(gdrive_location,"/dox_data/maxdrop.r2maf.filtered.csv.gz"), stringsAsFactors = F)
serie = serie %>% select(rsid, A1, A2, CHR, BP, BETA, P)

eqtl = read_qtls(paste0(DATADIR,"/panama_qq_boot_1e+06/"))  %>% left_join(snploc, by=c(cis_snp="snpid")) 

join_all = inner_join(eqtl, serie, by=c("RSID"="rsid")) 

#join_all = inner_join(eqtl %>% select(RSID), serie %>% select(rsid), by=c("RSID"="rsid")) 
chroms=sort(unique(eqtl$chr))
join_all=foreach(chrom=chroms, .combine = bind_rows) %do% { 
  print(chrom)
  eqtl_chr = eqtl %>% filter(chr==chrom)
  rownames(eqtl_chr) = eqtl_chr$RSID
  serie_chr = serie %>% filter(CHR==substr(chrom,4,10))
  rownames(serie_chr) = serie_chr$rsid
  shared_rsid=intersect(rownames(eqtl_chr),rownames(serie_chr))
  
}