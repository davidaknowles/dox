#!/usr/bin/env Rscript

# The loading order matters because they both have a function named `select`
suppressMessages(library("biomaRt"))
suppressMessages(library("dplyr"))

# Input files ------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
f_est_counts <- args[1]
f_tpm <- args[2]
files <- args[-1:-2]
# files <- c("/scratch/midway/jdblischak/burridge2016/con1-0um/abundance.tsv",
#            "/scratch/midway/jdblischak/burridge2016/con1-1um/abundance.tsv",
#            "/scratch/midway/jdblischak/burridge2016/con3-0um/abundance.tsv",
#            "/scratch/midway/jdblischak/burridge2016/con3-1um/abundance.tsv")
stopifnot(files > 0, file.exists(files))
num_files <- length(files)
# Extract sample names
fname_parts <- strsplit(dirname(files), "/")
sample_names <- sapply(fname_parts, function(x) x[length(x)])
# Input first file to learn the transcripts quantified
f1 <- read.delim(files[1], stringsAsFactors = FALSE)
num_transcripts <- nrow(f1)
annotation <- f1[, c("target_id", "length")]

# Create matrices to store results
mat_template <- matrix(nrow = num_transcripts, ncol = length(files),
                       dimnames = list(rownames(f1), sample_names))
mat_eff_length <- mat_template
mat_est_counts <- mat_template
mat_tpm <- mat_template

# Add first file
mat_eff_length[, 1] <- f1$eff_length
mat_est_counts[, 1] <- f1$est_counts
mat_tpm[, 1] <- f1$tpm

for (i in 2:num_files) {
  f <- read.delim(files[i], stringsAsFactors = FALSE)
  stopifnot(f[, c("target_id", "length")] == annotation)
  mat_eff_length[, i] <- f$eff_length
  mat_est_counts[, i] <- f$est_counts
  mat_tpm[, i] <- f$tpm
}

stopifnot(cor(mat_eff_length) > 0.99)

# Download gene info from Ensembl ----------------------------------------------

# The sequences available for download from kallisto website are from Ensembl
# release 79 in March 2015
ensembl <- useMart(host = "Mar2015.archive.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")

gene_info <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                                  "chromosome_name", "external_gene_name",
                                  "gene_biotype"),
                   mart = ensembl)

# The transcripts in the kallisto download are only a subset of the available
# biotypes:
# $ zcat transcriptome-ensembl-GRCh38.fa.gz | grep ">" | cut -d" " -f5 | sort | uniq -c | sed s/"gene_biotype:"//

biotypes_kallisto <- c("IG_C_gene",
                       "IG_C_pseudogene",
                       "IG_D_gene",
                       "IG_J_gene",
                       "IG_J_pseudogene",
                       "IG_V_gene",
                       "IG_V_pseudogene",
                       "polymorphic_pseudogene",
                       "processed_pseudogene",
                       "protein_coding",
                       "pseudogene",
                       "transcribed_processed_pseudogene",
                       "transcribed_unitary_pseudogene",
                       "transcribed_unprocessed_pseudogene",
                       "translated_processed_pseudogene",
                       "translated_unprocessed_pseudogene",
                       "TR_C_gene",
                       "TR_D_gene",
                       "TR_J_gene",
                       "TR_J_pseudogene",
                       "TR_V_gene",
                       "TR_V_pseudogene",
                       "unitary_pseudogene",
                       "unprocessed_pseudogene")

gene_info_kallisto <- gene_info[gene_info$gene_biotype %in% biotypes_kallisto, ]

# Sort everything
annotation <- annotation[order(annotation$target_id), ]
gene_info_kallisto <- gene_info_kallisto[order(gene_info_kallisto$ensembl_transcript_id), ]
stopifnot(annotation$target_id == gene_info_kallisto$ensembl_transcript_id)
mat_eff_length <- mat_eff_length[order(gene_info_kallisto$ensembl_transcript_id), ]
mat_est_counts <- mat_est_counts[order(gene_info_kallisto$ensembl_transcript_id), ]
mat_tpm <- mat_tpm[order(gene_info_kallisto$ensembl_transcript_id), ]

# Sum per transcript -----------------------------------------------------------

df_est_counts <- cbind(gene_info_kallisto, mat_est_counts)
df_est_counts_gene <- df_est_counts %>%
  select(-ensembl_transcript_id) %>%
  group_by(ensembl_gene_id, chromosome_name, external_gene_name, gene_biotype) %>%
  summarise_each(funs(sum)) %>%
  arrange(ensembl_gene_id)

df_tpm <- cbind(gene_info_kallisto, mat_tpm)
df_tpm_gene <- df_tpm %>%
  select(-ensembl_transcript_id) %>%
  group_by(ensembl_gene_id, chromosome_name, external_gene_name, gene_biotype) %>%
  summarise_each(funs(sum)) %>%
  arrange(ensembl_gene_id)

stopifnot(df_est_counts_gene$ensembl_gene_id == df_tpm_gene$ensembl_gene_id)

# Filter genes -----------------------------------------------------------------

# Only include protein_coding genes on main chromosomes
df_est_counts_gene <- df_est_counts_gene %>% filter(gene_biotype == "protein_coding",
                                                    chromosome_name %in% c(1:22, "X", "Y", "MT"))
df_tpm_gene <- df_tpm_gene %>% filter(gene_biotype == "protein_coding",
                                      chromosome_name %in% c(1:22, "X", "Y", "MT"))

# Output -----------------------------------------------------------------------

write.table(df_est_counts_gene, f_est_counts,
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(df_tpm_gene, f_tpm,
            quote = FALSE, sep = "\t", row.names = FALSE)
