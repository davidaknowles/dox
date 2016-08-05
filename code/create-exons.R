# Create exons file for for mapping reads to genes with Subread featureCounts.

# Usage:
#   Rscript create-exons.R > out.saf

# Notes:
# + This includes only coding genes, i.e. gene_biotype == "protein_coding"
# + Uses grch37 Ensembl archive
# + Output is in Simplified Annotation Format (SAF)
#     + Columns: GeneID, Chr, Start, End, Strand
#     + Coordinates are 1-based, inclusive on both ends
# + Contains duplicate and overlapping exons (featureCounts handles this)

# Download human exons (grch37, hg19)
library("biomaRt")
#ensembl <- useMart(host = "grch37.ensembl.org",
ensembl <- useMart(host = "dec2015.archive.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")
# attributePages(ensembl)
# [1] "feature_page" "structure" "homologs" "sequences" "snp" "snp_somatic"
atts <- listAttributes(ensembl, page = "feature_page")
# atts[grep("strand", atts$description, ignore.case = TRUE), ]
atts_struct <- listAttributes(ensembl, page = "structure")
# atts_struct[grep("exon", atts_struct$description, ignore.case = TRUE), ]
exons_all <- getBM(attributes = c("ensembl_gene_id", "ensembl_exon_id",
                                  "chromosome_name", "exon_chrom_start",
                                  "exon_chrom_end", "strand",
                                  "external_gene_name",
                                  "gene_biotype"),
                   mart = ensembl)
exons_final <- exons_all[exons_all$chromosome_name %in% c(1:22, "X", "Y", "MT") &
                         exons_all$gene_biotype == "protein_coding",
                         c("ensembl_gene_id", "chromosome_name", "exon_chrom_start",
                           "exon_chrom_end", "strand", "external_gene_name")]
colnames(exons_final) <- c("GeneID", "Chr", "Start", "End", "Strand", "Name")
# Sort by chromosome and position
exons_final <- exons_final[order(exons_final$Chr,
                                 exons_final$Start,
                                 exons_final$End), ]
# Fix chromosome names
exons_final$Chr <- paste0("chr", exons_final$Chr)
exons_final$Chr <- sub("chrMT", "chrM", exons_final$Chr)
# Fix strand
exons_final$Strand <- ifelse(exons_final$Strand == 1, "+", "-")

# Save as tab-separated file in Simplified Annotation Format (SAF)
write.table(exons_final, "", quote = FALSE, sep = "\t",
            row.names = FALSE)
