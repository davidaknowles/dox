#source("load_data.R")

pred_exp=foreach(i=1:22) %do% {
  cat(".")
  chrom=paste0("chr",i)
  fn=paste0(DATADIR,"predicted_expression_log/",chrom,".txt.gz")
  if (file.exists(fn)) { read.table(fn, sep="\t") %>% t() %>% as.data.frame() } else {NULL}
}

pred_exp=Reduce(rbind, pred_exp) # only way to keep rownames it seems

gz=gzfile(paste0(DATADIR,"predicted_log_cpm.txt.gz"),"w")
pred_exp %>% write.table(gz, quote=F, col.names=F, row.names = T, sep="\t")
close(gz)


#s=apply(input,1,sd)
#m=apply(input,1,mean)

#pred_exp = sweep(pred_exp, 1, s[rownames(pred_exp)], "*")
#pred_exp = sweep(pred_exp, 1, m[rownames(pred_exp)], "+")

pca <- prcomp(t(pred_exp), scale. = TRUE)
variances <- pca$sdev^2
explained <- variances / sum(variances)

pca_data <- cbind(anno, pca$x[, 1:5])

pca_data %>% mutate(conc=as.factor(conc)) %>% 
  ggplot( aes(x = PC1, y = PC2, group=individual, color=conc)) +
  geom_text(aes(label = conc)) +
  geom_path(aes(alpha=.3)) +
  labs(x = sprintf("PC%d (%.2f%%)", 1, round(explained[1] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 2, round(explained[2] * 100, 2))) # + xlim(75,100) + ylim(-25,25)


