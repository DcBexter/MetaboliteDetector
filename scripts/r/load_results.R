# Loads a results CSV exported from MetaboliteDetector (batch window -> "Export CSV...")
# as the data.frame `data`, like the old built-in R console.
# Usage: Rscript load_results.R results.csv   (or set `csv_path` and source this script)
args <- commandArgs(trailingOnly = TRUE)
if (!exists("csv_path")) csv_path <- if (length(args)) args[1] else "results.csv"

# Columns: Substance, RT_s, Quant_ions, then one column per run (peak area sum; NA = not found).
# The file starts with a byte order mark (for Excel); text that a spreadsheet would read as a
# formula (starting with = + - @) has a leading ' that is taken off again here.
unguard <- function(x) sub("^'(?=[=+@\t\r-])", "", x, perl = TRUE)
data <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
                 colClasses = c(Substance = "character", Quant_ions = "character"))
names(data) <- unguard(names(data))
data$Substance <- unguard(data$Substance)
data[-(1:3)] <- lapply(data[-(1:3)], as.numeric)  # an all-NA run column would read as logical
stopifnot(all(c("Substance", "RT_s", "Quant_ions") %in% names(data)))
runs <- data[, -(1:3), drop = FALSE]

# Replicate groups of the runs (exported with --replicates or a regex in the batch window):
# `groups` is a named character vector, run name -> group, in the order of the run columns.
groups_csv <- paste0(tools::file_path_sans_ext(csv_path), ".groups.csv")
groups <- NULL
if (file.exists(groups_csv)) {
  g <- read.csv(groups_csv, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
               colClasses = "character", na.strings = character(0))
  groups <- setNames(unguard(g$Group), unguard(g$Run))[names(runs)]
}

# Settings used (`Setting,Value`, written next to every export) as the named character vector
# `settings`, and the integral of every quant ion (exported with "Export quant ions" or
# md-cli --ion-integrals; Substance, Quant_ion, then one column per run) as the data.frame `ions`.
settings_csv <- paste0(tools::file_path_sans_ext(csv_path), ".settings.csv")
settings <- NULL
if (file.exists(settings_csv)) {
  s <- read.csv(settings_csv, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
                colClasses = "character", na.strings = character(0))
  settings <- setNames(unguard(s$Value), s$Setting)
}
ions_csv <- paste0(tools::file_path_sans_ext(csv_path), ".ions.csv")
ions <- NULL
if (file.exists(ions_csv)) {
  ions <- read.csv(ions_csv, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
                   colClasses = c(Substance = "character"))
  names(ions) <- unguard(names(ions))
  ions$Substance <- unguard(ions$Substance)
  ions[-(1:2)] <- lapply(ions[-(1:2)], as.numeric)
}

str(data)
