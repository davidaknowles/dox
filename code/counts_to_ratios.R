ratios = counts %>% 
  mutate(clu = str_split_fixed(rownames(counts), ":", 4)[,4]) %>%
  group_by(clu) %>% 
  mutate_all( funs( ./sum(.) ) ) %>% 
  ungroup() %>%
  as.data.frame() %>% 
  set_rownames(rownames(counts)) %>% 
  select(-clu)

ratios = ratios[rowMeans(is.na(ratios)) <= 0.4,,drop=F ]
row_means = rowMeans(ratios, na.rm = T)
row_means_outer = outer(row_means, rep(1,ncol(ratios)))
ratios[is.na(ratios)] = row_means_outer[is.na(ratios)]