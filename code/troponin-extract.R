# /usr/bin/env Rscript

# The cardiac tropnonin measurements are in 3 Excel files. This script extracts
# the data into a tab-separated, plain-text file. Execute the script from the
# code subdirectory.

# In each Excel file, each sheet represents one 96-well plate. The plate is read
# on the spectrometer twice. Thus each sheet has 3 96-well configurations: the
# values from the first reading, the values from the second reading, and the
# identify of the wells.

# cTnI_ELISA_8-12-15.xlsx - 2 plates
# cTnI_ELISA-10-26-15.xlsx - 2 plates
# cTnI_ELISA_11-5-15.xlsx - 3 plates

# Note: I (John) had to make a few small hand edits to the names of the wells so
# that they were consistent across the experiments.

# Setup ------------------------------------------------------------------------

library("openxlsx")
library("reshape2")
library("plyr")
library("stringr")

data_dir <- "../data"

xlsx_files <- file.path(data_dir, c("cTnI_ELISA_8-12-15.xlsx",
                                    "cTnI_ELISA-10-26-15.xlsx",
                                    "cTnI_ELISA_11-5-15.xlsx"))
stopifnot(file.exists(xlsx_files))

extract_plate <- function(x) {
  # Extract the absorbance values from one Excel sheet
  #
  # x - data frame returned by openxlsx::read.xlsx
  #
  stopifnot(is.data.frame(x), nrow(x) >= 26, ncol(x) == 14)
  colnames(x)[1] <- "row"
  # Remove column that specifies the absorbance wavelength
  x <- x[, -ncol(x)]
  # Isolate values of first reading
  reading1 <- x[1:8, ]
  reading1_long <- melt(reading1, id.vars = "row",
                        variable.name = "column", value.name = "absorbance")
  # Isolate values of second reading
  reading2 <- x[10:17, ]
  reading2_long <- melt(reading2, id.vars = "row",
                        variable.name = "column", value.name = "absorbance")
  # Isolate identities of each well
  identities <- x[19:26, ]
  identities$row <- LETTERS[1:8]
  identities_long <- melt(identities, id.vars = "row",
                          variable.name = "column", value.name = "name")
  # Empty wells were assigned the identity "x"
  identities_long$name[identities_long$name == "x"] <- "empty"
  stopifnot(nrow(reading1_long) == 96,
            nrow(reading2_long) == 96,
            nrow(identities_long) == 96,
            reading1_long$row == rep(LETTERS[1:8], 12),
            reading1_long$column == rep(1:12, each = 8),
            reading2_long$row == rep(LETTERS[1:8], 12),
            reading2_long$column == rep(1:12, each = 8),
            identities_long$row == rep(LETTERS[1:8], 12),
            identities_long$column == rep(1:12, each = 8))
  # Combine the 3 sets of data
  reading1_long$reading <- 1
  reading2_long$reading <- 2
  combined <- rbind(reading1_long, reading2_long)
  combined$name <- identities_long$name
  return(combined)
}

# Extract data from xlsx files -------------------------------------------------

# Iterate through each sheet of each file to extract absorbance values. Store
# each plate in a list.
final <- list()
for (f in xlsx_files) {
  for (sheet in 1:3) {
  raw <- read.xlsx(f, sheet = sheet)
  if (is.null(raw)) next
  plate_name <- paste(basename(f), "sheet", sheet, sep = ".")
  final[[plate_name]] <- extract_plate(raw)
  }
}

# Format data ------------------------------------------------------------------

# Combine all the plates into one data frame
final_df <- ldply(final, .id = "experiment")
stopifnot(!anyNA(final_df), nrow(final_df) %% 96 == 0)

# Format names
# Most of the standards are named with ng/mL, but some have just ng.
for (i in 1:nrow(final_df)) {
  name_len <- nchar(final_df$name[i])
  if (str_sub(final_df$name[i], start = name_len - 1, end = name_len) == "ng") {
    final_df$name[i] <- paste0(final_df$name[i], "/mL")
  }
}
# Add a column type, which specifies the type of the well
final_df$type <- ifelse(grepl("ng/mL", final_df$name), "standard", "unknown")
final_df$type[final_df$name == "water"] <- "control"
final_df$type[final_df$name == "galactose"] <- "control"
final_df$type[final_df$name == "empty"] <- "control"
# table(final_df$type)
# Extract the cell_line and dosage from the names of the unknown wells
name_split <- str_split_fixed(final_df$name, pattern = " ", n = 2)
name_split[final_df$type != "unknown", ] <- NA
final_df$cell_line <- name_split[, 1]
final_df$cell_line <- str_replace(final_df$cell_line, "-", ".")
final_df$cell_line[!is.na(final_df$cell_line)] <- paste0("c",
                          final_df$cell_line[!is.na(final_df$cell_line)])
final_df$dosage <- name_split[, 2]
final_df$dosage[final_df$dosage == "control"] <- "C"
# cTnI_ELISA_11-5-15.xlsx used "c" instead of "control"
final_df$dosage[final_df$dosage == "c"] <- "C"
# Extract the concentrations of the standards
final_df$known_conc <- ifelse(final_df$type == "standard",
                              str_replace(final_df$name, "ng/mL", ""),
                              NA)

# Save data --------------------------------------------------------------------

write.table(final_df, file = file.path(data_dir, "troponin.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
