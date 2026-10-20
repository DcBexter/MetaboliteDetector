# Loads a results CSV exported from MetaboliteDetector (batch window -> "Export CSV...")
# as the data.frame `data`, like the old built-in R console.
# Usage: Rscript load_results.R results.csv   (or set `csv_path` and source this script)
args <- commandArgs(trailingOnly = TRUE)
if (!exists("csv_path")) csv_path <- if (length(args)) args[1] else "results.csv"

# Columns: Substance, RT_s, Quant_ions, then one column per run (peak area sum; NA = not found).
data <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8",
                 colClasses = c(Substance = "character", Quant_ions = "character"))
data[-(1:3)] <- lapply(data[-(1:3)], as.numeric)  # an all-NA run column would read as logical
stopifnot(all(c("Substance", "RT_s", "Quant_ions") %in% names(data)))
runs <- data[, -(1:3), drop = FALSE]

str(data)
