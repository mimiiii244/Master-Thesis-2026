setwd("/Users/amyhu/Research/AI_Energy_Thesis_Code")
# 04 Generation vs Utilities
# Change the working-directory line above for another computer.
data_dir <- "09_Shared_Data/inputs"
output_dir <- "04_Generation_vs_Utilities/outputs"
for (folder_name in c("tables", "models", "figures", "data")) {
  dir.create(file.path(output_dir, folder_name), recursive = TRUE, showWarnings = FALSE)
}
source("10_Shared_Functions/Shared_Functions.R")
library(BigVAR)

# Generation versus regulated utilities
# Read the included verified matched-network input. The 29-stock Technology and
# 13-stock hardware baskets are sender columns in the same 90-stock system.
# 1,596 synchronized fixed-penalty bootstrap refits.

period_models <- readRDS("09_Shared_Data/frozen_models/matched_network_models.rds")
matched_tickers <- colnames(period_models$post$y)

technology_groups <- list(DY_Technology_29 = intersect(sector_groups$Technology, matched_tickers),
  Original_CHIP_hardware_13 = intersect(chip_tickers, matched_tickers))

electricity_groups <- list(Generation_3 = intersect(sector_groups$Generation, matched_tickers),
  Regulated_Utilities_19 = intersect(utility_tickers, matched_tickers), Merchant_VST_NRG = c("VST",
    "NRG"), EXC_separate = "EXC")

electricity_firms <- c(electricity_groups[[1]], electricity_groups[[2]])

stopifnot(length(technology_groups[[1]]) == 29, length(technology_groups[[2]]) == 13,
  length(electricity_firms) == 22, all(technology_groups[[2]] %in% technology_groups[[1]]))

definition_lists <- list(senders = technology_groups, electric = electricity_groups)
for (kind in names(definition_lists)) {
  a <- definition_lists[[kind]]
  write.csv(do.call(rbind, lapply(names(a), function(g) {
    data.frame(definition = g, ticker = a[[g]])
  })), file.path(output_dir, "tables", paste0(kind, "_definitions", ".csv")), row.names = FALSE)
}

measure_electricity_dependence <- function(o) {
  y <- o$y
  C <- cor(y)
  out <- list()
  add <- function(level, entity, source, metric, value) {
    out[[length(out) + 1]] <<- data.frame(level, entity, source, metric, value)
  }
  for (src in names(technology_groups)) {
    jj <- technology_groups[[src]]
    tp <- rowMeans(y[, jj, drop = FALSE])
    fc <- cor(y[, electricity_firms, drop = FALSE], tp)[, 1]
    pair <- rowMeans(C[electricity_firms, jj, drop = FALSE])
    alloc <- 100 * rowSums(o$D[electricity_firms, jj, drop = FALSE])
    zz <- 100 * rowSums(o$Z[electricity_firms, jj, drop = FALSE])
    vals <- list(DY_technology_to_firm_pct = alloc, zero_technology_to_firm_pct = zz,
      signed_portfolio_correlation = fc, mean_pairwise_signed_correlation = pair)
    for (metric in names(vals)) {
      for (f in electricity_firms) {
        add("firm", f, src, metric, vals[[metric]][f])
      }
    }
    for (g in names(electricity_groups)) {
      ii <- electricity_groups[[g]]
      for (metric in names(vals)) {
        add("group", g, src, metric, mean(vals[[metric]][ii]))
      }
      add("group", g, src, "equal_weight_group_portfolio_correlation", cor(rowMeans(y[,
        ii, drop = FALSE]), tp))
      add("group", g, src, "DY_electricity_to_technology_pct", 100 * mean(rowSums(o$D[jj,
        ii, drop = FALSE])))
      add("group", g, src, "DY_electricity_to_technology_per_sender_pct", 100 *
        mean(o$D[jj, ii, drop = FALSE]))
    }
  }
  a <- do.call(rbind, out)
  a$key <- paste(a$level, a$entity, a$source, a$metric, sep = "|")
  a
}

points <- do.call(rbind, lapply(names(period_models), function(tag) {
  x <- measure_electricity_dependence(period_models[[tag]])
  x$period <- tag
  x
}))

write.csv(points, file.path(output_dir, "tables", "point_estimates.csv"), row.names = FALSE)

main <- points[points$period == "pre_main", ]
post <- points[points$period == "post", ]

stopifnot(identical(main$key, post$key))

change <- main[, c("level", "entity", "source", "metric", "key")]
change$pre <- main$value
change$post <- post$value
change$delta <- post$value - main$value

write.csv(change, file.path(output_dir, "tables", "main_changes.csv"), row.names = FALSE)

generation_utility_contrast <- function(x) {
  a <- x[x$level == "group" & x$entity == "Generation_3", ]
  b <- x[x$level == "group" & x$entity == "Regulated_Utilities_19", ]
  stopifnot(identical(a$source, b$source), identical(a$metric, b$metric))
  data.frame(source = a$source, metric = a$metric, value = a$value - b$value, key = paste(a$source,
    a$metric, sep = "|"))
}

cc <- do.call(rbind, lapply(names(period_models), function(tag) {
  x <- generation_utility_contrast(points[points$period == tag, ])
  x$period <- tag
  x
}))

write.csv(cc, file.path(output_dir, "tables", "generation_minus_utilities_levels.csv"),
  row.names = FALSE)

pre <- as.matrix(read.csv(file.path(data_dir, "R_pre_ai_returns.csv"), row.names = 1,
  check.names = FALSE))
pr <- as.matrix(read.csv(file.path(data_dir, "R_ai_energy_lasso_100_daily_returns.csv"),
  row.names = 1, check.names = FALSE))
ff <- read.csv(file.path(data_dir, "FF3_daily_R.csv"))

data <- list(pre_main = list(r = pre[, matched_tickers], f = ff[match(rownames(pre), ff$date),
  ]), post = list(r = pr[, matched_tickers], f = ff[match(rownames(pr), ff$date), ]))

ident <- list()

for (tag in names(data)) {
  for (src in names(technology_groups)) {
    d <- data[[tag]]
    a <- adjust_ff_factors(matrix(rowMeans(d$r[, technology_groups[[src]]]), ncol = 1),
      d$f)
    b <- rowMeans(period_models[[tag]]$y[, technology_groups[[src]]])
    gap <- max(abs(a - b))
    stopifnot(gap < 1e-12)
    ident[[length(ident) + 1]] <- data.frame(period = tag, source = src, max_gap = gap)
  }
}

write.csv(do.call(rbind, ident), file.path(output_dir, "tables", "portfolio_residual_identity.csv"),
  row.names = FALSE)

B <- 399L

jobs <- expand.grid(period = c("pre_main", "post"), block = c(10, 20), stringsAsFactors = FALSE)

run_generation_bootstrap <- function(j) {
  tag <- jobs$period[j]
  L <- jobs$block[j]
  d <- data[[tag]]
  set.seed(20261016 + 1000 * j + 71)
  values <- matrix(NA_real_, B, nrow(main), dimnames = list(NULL, main$key))
  contrasts <- matrix(NA_real_, B, nrow(generation_utility_contrast(main)), dimnames = list(NULL,
    generation_utility_contrast(main)$key))
  for (b in seq_len(B)) {
    ix <- moving_block_indices(nrow(d$r), L)
    y <- adjust_ff_factors(d$r[ix, , drop = FALSE], d$f[ix, , drop = FALSE])
    o <- fit_fixed_penalty(y, period_models[[tag]]$lambda)
    a <- measure_electricity_dependence(o)
    stopifnot(identical(a$key, main$key))
    values[b, ] <- a$value
    contrasts[b, ] <- generation_utility_contrast(a)$value
  }
  saveRDS(list(values = values, contrasts = contrasts, period = tag, block = L), paste0(file.path(output_dir,
    "bootstrap_"), tag, "_L", L, ".rds"))
  list(values = values, contrasts = contrasts)
}

# Run in the current R session: forked workers can disappear in macOS RStudio.
# Each job sets its original seed, so sequential execution preserves its draws.
boots <- lapply(seq_len(nrow(jobs)), function(j) {
  message("Generation bootstrap: ", jobs$period[j], ", block ", jobs$block[j])
  result <- run_generation_bootstrap(j)
  stopifnot(is.list(result), is.matrix(result$values),
    is.matrix(result$contrasts),
    identical(dim(result$values), c(B, nrow(main))),
    identical(colnames(result$values), main$key),
    identical(dim(result$contrasts), c(B, nrow(generation_utility_contrast(main)))),
    identical(colnames(result$contrasts), generation_utility_contrast(main)$key),
    all(is.finite(result$values)), all(is.finite(result$contrasts)))
  result
})

changes <- list()
did <- list()
cpre <- generation_utility_contrast(main)
cpost <- generation_utility_contrast(post)

for (L in c(10, 20)) {
  i <- which(jobs$period == "pre_main" & jobs$block == L)
  j <- which(jobs$period == "post" & jobs$block == L)
  b <- boots[[j]]$values - boots[[i]]$values
  v <- change
  v$block <- L
  v$low <- apply(b, 2, quantile, 0.025)
  v$high <- apply(b, 2, quantile, 0.975)
  v$B <- B
  changes[[as.character(L)]] <- v
  b <- boots[[j]]$contrasts - boots[[i]]$contrasts
  v <- cpre[, c("source", "metric", "key")]
  v$pre_gap <- cpre$value
  v$post_gap <- cpost$value
  v$contrast_change <- cpost$value - cpre$value
  v$low <- apply(b, 2, quantile, 0.025)
  v$high <- apply(b, 2, quantile, 0.975)
  v$block <- L
  v$B <- B
  did[[as.character(L)]] <- v
}

write.csv(do.call(rbind, changes), file.path(output_dir, "tables", "firm_and_group_change_intervals.csv"),
  row.names = FALSE)

write.csv(do.call(rbind, did), file.path(output_dir, "tables", "generation_vs_utilities_change_contrast.csv"),
  row.names = FALSE)

write.csv(data.frame(period = jobs$period, block = jobs$block, B = vapply(boots, function(x) {
  nrow(x$values)
}, integer(1))), file.path(output_dir, "tables", "bootstrap_completion.csv"), row.names = FALSE)

# Firm tables, solver reconstruction and figures
# Differentiate group-index correlations from means of firm-level correlations.
# Keep the original 0.001 percentage-point direct-solver comparison bound.

x <- read.csv(file.path(output_dir, "tables/main_changes.csv"))

ci <- read.csv(file.path(output_dir, "tables/generation_vs_utilities_change_contrast.csv"))

points <- read.csv(file.path(output_dir, "tables/point_estimates.csv"))
fci <- read.csv(file.path(output_dir, "tables/firm_and_group_change_intervals.csv"))

for (src in unique(x$source)) {
  a <- x[x$level == "firm" & x$source == src & x$metric == "DY_technology_to_firm_pct",
    ]
  b <- x[x$level == "firm" & x$source == src & x$metric == "signed_portfolio_correlation",
    ]
  b <- b[match(a$entity, b$entity), ]
  t <- data.frame(firm = a$entity, group = ifelse(a$entity %in% c("VST", "NRG"), "Merchant_pair",
    ifelse(a$entity == "EXC", "EXC_separate", "Regulated_Utilities")), DY_pre = a$pre,
    DY_post = a$post, DY_change = a$delta, Corr_pre = b$pre, Corr_post = b$post, Corr_change = b$delta)
  write.csv(t, file.path(output_dir, "tables", paste0("firm_comparison_", src, ".csv")),
    row.names = FALSE)
}

grDevices::cairo_pdf(file.path(output_dir, "figures/firm_changes.pdf"), width = 12, height = 9,
  onefile = TRUE, family = "DejaVu Sans")

for (metric in c("DY_technology_to_firm_pct", "signed_portfolio_correlation")) {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  for (src in unique(x$source)) {
    a <- fci[fci$level == "firm" & fci$source == src & fci$metric == metric & fci$block ==
      20, ]
    yy <- rev(seq_len(nrow(a)))
    cols <- ifelse(a$entity %in% c("VST", "NRG", "EXC"), "#ba583d", "#2f668d")
    plot(a$delta, yy, xlim = range(a$low, a$high, 0), ylim = c(0.5, nrow(a) + 0.5),
      yaxt = "n", pch = 19, col = cols, xlab = if (metric == "DY_technology_to_firm_pct") {
        "Post-minus-pre DY allocation (pp)"
      }
      else {
        "Post-minus-pre portfolio correlation"
      }, ylab = "", main = if (src == "DY_Technology_29") {
        "Broad DY Technology (29)"
      }
      else {
        "Original chip/hardware (13)"
      })
    segments(a$low, yy, a$high, yy, col = cols)
    axis(2, yy, a$entity, las = 1, cex.axis = 0.8)
    abline(v = 0, lty = 2)
    abline(h = 19.5, lty = 3, col = "grey")
    mtext("95% conditional intervals; L=20; 399 draws", side = 3, line = 0.2, cex = 0.8)
  }
}

dev.off()

