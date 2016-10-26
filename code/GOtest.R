require(GO.db)
require(GOstats)
require(biomaRt)

GOtest <- function(my_ensembl_gene_universe,my_ensembl_gene_test,
                   pval_cutoff = .01,ontology=c("BP","CC","MF"),
                   conditional.method = F) {

  #ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
                     host="www.ensembl.org")
  
  # Convert Ensembl IDs to Entrez IDs
  ensembl_to_entrez <- getBM(attributes=c('ensembl_gene_id', 'entrezgene'), 
                             filters = 'ensembl_gene_id', 
                             values = my_ensembl_gene_universe, mart = ensembl)
  ii_entrez_gene <- which(ensembl_to_entrez$ensembl_gene_id %in% my_ensembl_gene_test)
  my_entrez_gene <- ensembl_to_entrez[ii_entrez_gene, ]
  
  ensembl_to_entrez <- unique(ensembl_to_entrez)
  my_entrez_gene <- unique(my_entrez_gene)

  hg_over <- setNames( foreach(ont=ontology) %do% {
    params <- new("GOHyperGParams",
                  geneIds=my_entrez_gene$entrezgene,
                  universeGeneIds = ensembl_to_entrez$entrezgene,
                  ontology=ontology,
                  pvalueCutoff=pval_cutoff,
                  conditional = conditional.method,
                  testDirection="over",
                  annotation="org.Hs.eg.db")
    hyperGTest(params)
  }, ontology)
  list(GO = hg_over, geneIds = ensembl_to_entrez )
}


plotHeatmap <- function(heatmapData, mar =  c(0, 0, 0, 0),
                        oma = c(0, 0, 0, 0), keysize = .5,
                        labCol = "" ) {
  require(broman)
  crayon <- brocolors("crayon")
  cols =  colorRampPalette(c("white", crayon["Banana Mania"],
                             crayon["Laser Lemon"],
                             crayon["Burnt Orange"], crayon["Orange Red"],
                             crayon["Plum"]))
  par( mfrow = c(1,1), mar = mar, oma = oma )
  require(gplots)
  mat <- matrix(heatmapData[ , -c(1:2)], ncol = 1)
  sepwidth <- c(.1, .01)
  sepcol <- "grey30"
  
  if (NCOL(mat) == 1) {
    mat <- cbind(mat, mat)
    sepwidth <- c(0, .01)
    sepcol <- "white"
  }
  rownames(mat) <- heatmapData$Term
  hmap <- heatmap.2(as.matrix(-log10(mat)), main = "",
                    labCol = labCol, cexCol = 1.2,
                    cexRow = 1.3,
                    sepcol = "white", sepwidth=sepwidth,
                    colsep = 1:ncol(mat), Colv = FALSE,
                    rowsep = 1:nrow(mat), dendrogram = "row",
                    srtCol = 0, adjCol = c(.5,NA), offsetCol = 1,
                    keysize = keysize, lwid = c(1.5, 7), lhei = c(.5, 5),
                    col = cols, trace ="none", key.title = "-log p-value" )
  return(hmap)
}


goHeatmapData <- function(golist){
  nlist <- length(golist)
  goIDlist <- lapply( 1:nlist, function(i) {
    golist[[i]][ , 1]
  })
  names(goIDlist) <- names(golist)
  
  goTermlist <- lapply( 1:nlist, function(i) {
    golist[[i]]$Term
  })
  
  goset <- unique( unlist(goIDlist) )
  gosetTerm <- unique( unlist(goTermlist) )
  
  heatmap_allgos <- data.frame(ID = goset)
  heatmap_allgos$Term <- gosetTerm
  
  pvals <- lapply(1:nlist, function(i) {
    foo <- rep(1-10^(-6), length(gosetTerm))
    foo[ which(goset %in% golist[[i]][ , 1]) ] <- golist[[i]]$Pvalue[ which( golist[[i]][ , 1] %in% goset )]
    return(foo)
  })
  pvals <- do.call(cbind, pvals)
  colnames(pvals) <- names(golist)
  heatmap_allgos <- cbind(heatmap_allgos, pvals)
  colnames(heatmap_allgos)[1] <- colnames(golist[[1]])[1]
  
  ind <- which( is.na(heatmap_allgos), arr.ind = TRUE)
  heatmap_allgos[ ind ] <- 1
  rownames(heatmap_allgos) <- heatmap_allgos[ , 1]
  
  iisig <- pvals > .01
  iisig_gos <- which(rowSums(iisig)!=nlist)
  
  heatmap_allgos <- heatmap_allgos[iisig_gos, ]
  
  return(heatmap_allgos)
}
