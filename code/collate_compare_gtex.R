require(data.table)
require(dplyr)
require(foreach)

pi0_files=dir()

foreach(f=pi0_files, .combine = bind_rows) %do% {
  read.table(f, header=T, stringsAsFactors = F)
} %>% write.table("gtex_pi0.txt", sep="\t", row.names=F, quote=F)