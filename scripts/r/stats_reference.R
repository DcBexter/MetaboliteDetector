# Reference values for crates/md-core/src/stats.rs tests (issue #24), computed with base R.
# Usage: Rscript stats_reference.R   (prints the Rust literals pasted into the tests)
# Table: 6 substances x 10 runs: replicate groups c1, c2 (runs 1-4, "ctrl"), d1, d2 (runs 5-8, "case") and the pool P.
set.seed(24)
grp <- c("ctrl", "ctrl", "ctrl", "ctrl", "case", "case", "case", "case", "P", "P")
D <- matrix(round(rlnorm(6 * 10, meanlog = 5, sdlog = 0.6), 1), nrow = 6)
D[2, 3] <- NA  # not found; the old table stores 0
D[5, ] <- D[5, 1]  # a constant substance: no spread, dropped by the PCA, p = 1
fmt <- function(x) paste0("[", paste(sprintf("%.17g", x), collapse = ", "), "]")
cat("data:\n"); for (i in 1:6) cat(" ", paste(ifelse(is.na(D[i, ]), "None", sprintf("Some(%s)", D[i, ])), collapse = ", "), "\n")
Z <- D; Z[is.na(Z)] <- 0  # missing as zero
a <- which(grp == "ctrl"); b <- which(grp == "case")
cat("welch independent (missing as 0):", fmt(sapply(1:6, function(i) if (sd(Z[i, a]) == 0 || sd(Z[i, b]) == 0) 1 else t.test(Z[i, a], Z[i, b])$p.value)), "\n")
cat("welch independent (missing skipped):", fmt(sapply(1:6, function(i) if (sd(D[i, a], na.rm = TRUE) == 0 || sd(D[i, b]) == 0) 1 else t.test(D[i, a], D[i, b])$p.value)), "\n")
# group means: one value per replicate group is only possible with several groups, so split the
# runs of each side into replicate groups of two runs: g1 = (1,2) (3,4); g2 = (5,6) (7,8)
m <- function(i, idx) mean(Z[i, idx])
gm1 <- function(i) c(m(i, 1:2), m(i, 3:4)); gm2 <- function(i) c(m(i, 5:6), m(i, 7:8))
cat("welch group means:", fmt(sapply(1:6, function(i) if (sd(gm1(i)) == 0 || sd(gm2(i)) == 0) 1 else t.test(gm1(i), gm2(i))$p.value)), "\n")
cat("paired:", fmt(sapply(1:6, function(i) tryCatch(t.test(Z[i, a], Z[i, b], paired = TRUE)$p.value, error = function(e) NaN))), "\n")
cat("means ctrl/case:", fmt(rowMeans(Z[, a])), fmt(rowMeans(Z[, b])), "\n")
cat("sd ctrl/case:", fmt(apply(Z[, a], 1, sd)), fmt(apply(Z[, b], 1, sd)), "\n")
cat("five numbers row 1 ctrl:", fmt(c(min(Z[1, a]), quantile(Z[1, a], c(.25, .5, .75), type = 7), max(Z[1, a]))), "\n")
# PCA like the old showPCA: drop the pool runs, z-score each substance over the runs, covariance eigen
s <- which(grp != "P")
keep <- apply(Z[, s], 1, sd) > 0
Zs <- t(apply(Z[keep, s], 1, function(v) (v - mean(v)) / sd(v)))
pc <- prcomp(t(Zs), center = TRUE, scale. = FALSE)
cat("kept substances:", which(keep), "\n")
cat("proportion:", fmt(pc$sdev^2 / sum(pc$sdev^2)), "\n")
cat("abs scores PC1:", fmt(abs(pc$x[, 1])), "\n")
cat("abs scores PC2:", fmt(abs(pc$x[, 2])), "\n")
cat("abs loadings PC1:", fmt(abs(pc$rotation[, 1])), "\n")
cat("abs loadings PC2:", fmt(abs(pc$rotation[, 2])), "\n")
