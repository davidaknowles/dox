#!/usr/bin/Rscript

# Check concordance of genotype calls with reference genome. Results written to
# standard error.
#
# Usage:
#
# Rscript convert-cgi-to-vcf-test.R vcf
#
# vcf - path to a VCF file
#
# Also see the preceding script convert-cgi-to-vcf.py.

suppressPackageStartupMessages(library("BSgenome"))
suppressPackageStartupMessages(library("SNPlocs.Hsapiens.dbSNP144.GRCh38"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("VariantAnnotation"))

# Helpful Bioconductor tutorial:
# http://bioconductor.org/help/course-materials/2014/BioC2014/Lawrence_Tutorial.pdf

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
  stop("Incorrect number of arguments.",
       "\n\nUsage: Rscript convert-cgi-to-vcf-test.R vcf")
vcf_fname <- args[1]
# vcf_fname <- "/mnt/gluster/home/jdblischak/ober/vcf/dox-hg38-chr21.vcf.gz"
stopifnot(file.exists(vcf_fname))
vcf <- readVcf(vcf_fname, genome = "hg38")
# Checking distribution of number of samples with genotype information per locus
# ns <- info(vcf)$NS
# summary(ns)
# hist(ns)

# Logging
log_all <- 0            # All SNPs in VCF file
log_rsid <- 0           # SNPs with an rsID
log_unavailable <- 0    # SNPs not found in reference genome
log_loc <- 0            # SNPs with mismatched location
log_allele <- 0         # SNPs with mismatched alleles

snps_vcf <- granges(vcf)
stopifnot(class(snps_vcf) == "GRanges")
log_all <- length(snps_vcf)

# Subset to those SNPs which have an rsID
snps_vcf_rs <- snps_vcf[str_sub(names(snps_vcf), 1, 2) == "rs"]
stopifnot(length(snps_vcf_rs) > 0,
          length(snps_vcf_rs) < length(snps_vcf))
log_rsid <- length(snps_vcf_rs)

# Unforunately there is no way to handle missing SNPs, so I do this 1-by-1
for (rsid in names(snps_vcf_rs)) {
  message("Processing:\t", rsid)
  gr_vcf <- snps_vcf_rs[rsid]
  gr_ref <- tryCatch(snpid2grange(SNPlocs.Hsapiens.dbSNP144.GRCh38,
                                  snpid = rsid),
                     error = function(e) return(FALSE))
  if (class(gr_ref) == "GRanges") {
    stopifnot(length(gr_ref) == 1)
    start_vcf <- start(snps_vcf_rs[rsid])
    start_ref <- start(gr_ref)
    # Do the locations match?
    if (start_vcf != start_ref) {
      message("Inconsistent location:\t", rsid, "\tvcf:\t", start_vcf,
              "\tref:\t", start_ref)
      log_loc <- log_loc + 1
    }
    ref_vcf <- as.character(mcols(gr_vcf)$REF)
    alt_vcf <- as.character(mcols(gr_vcf)$ALT[[1]])
    # They only return the ambiguity code for the reference genotype. So
    # annoying.
    alleles_vcf <- paste(sort(c(ref_vcf, alt_vcf)), collapse = "")
    ambiguity_ref <- mcols(gr_ref)$alleles_as_ambig
    alleles_ref <- IUPAC_CODE_MAP[ambiguity_ref]
    # Do the alleles match?
    if (alleles_vcf != alleles_ref) {
      message("Inconsistent alleles:\t", rsid, "\tvcf:\t", alleles_vcf,
              "\tref:\t", alleles_ref)
      log_allele <- log_allele + 1
    }
  } else {
    message("Unavailable:\t", rsid)
    log_unavailable <- log_unavailable + 1
  }
}

message("\nLog:")
message("All SNPs in VCF file:\t", log_all)
message("SNPs with an rsID:\t", log_rsid)
message("SNPs not found in reference genome:\t", log_unavailable)
message("SNPs with mismatched location:\t", log_loc)
message("SNPs with mismatched alleles:\t", log_allele)
