setwd("/Users/amyhu/Research/AI_Energy_Thesis_Code")
# 05 Portfolio Analysis
# Change the working-directory line above for another computer.
data_dir <- "09_Shared_Data/inputs"
output_dir <- "05_Portfolio_Analysis/outputs"
for (folder_name in c("tables", "models", "figures", "data")) {
  dir.create(file.path(output_dir, folder_name), recursive = TRUE, showWarnings = FALSE)
}
source("10_Shared_Functions/Shared_Functions.R")
library(sandwich)
library(lmtest)

# Portfolio exposures and bootstrap inference
# Equal-weight log-return indices; original complete-case control alignment.
# 1,000/1,200 draws, all seeds, winsorisation and normal-tail HAC conventions.

pre_returns <- read.csv(file.path(data_dir, "R_pre_ai_returns.csv"), row.names = 1, check.names = FALSE)
post_returns <- read.csv(file.path(data_dir, "R_post_ai_returns.csv"), row.names = 1,
  check.names = FALSE)
common <- sort(intersect(colnames(pre_returns), colnames(post_returns)))
pre_returns <- pre_returns[, common]
post_returns <- post_returns[, common]

stopifnot(length(common) == 90)

market_controls <- readRDS("09_Shared_Data/frozen_models/market_controls.rds")
align_control_dates <- function(r) {
  d <- data.frame(Date = as.Date(rownames(r)), r, check.names = FALSE)
  d <- merge(d, market_controls$sp, by = "Date")
  d <- merge(d, market_controls$dg, by = "Date", all.x = TRUE)
  d <- d[complete.cases(d), ]
  d[order(d$Date), ]
}

construct_portfolios <- function(r) {
  d <- align_control_dates(r)
  d$oil <- rowMeans(d[, c("XOM", "CVX")])
  d$chip <- rowMeans(d[, chip_tickers])
  d$hyper <- rowMeans(d[, hyperscaler_tickers])
  d$merchant <- rowMeans(d[, c("VST", "NRG")])
  d$utility <- rowMeans(d[, utility_tickers])
  d
}

P0 <- construct_portfolios(pre_returns)
P1 <- construct_portfolios(post_returns)

pcnames <- c("merchant-chip", "VST-chip", "NRG-chip", "utility-chip", "merchant-hyper",
  "utility-hyper", "diff_merch_util_chip")

portfolio_statistics <- function(d, trim = FALSE) {
  q <- as.matrix(d[, c("mkt", "oil", "dy10", "chip", "hyper", "merchant", "VST", "NRG",
    "utility")])
  if (trim) {
    for (j in c(1, 2, 4:9)) {
      b <- quantile(q[, j], c(0.01, 0.99))
      q[, j] <- pmin(b[2], pmax(b[1], q[, j]))
    }
  }
  e <- factor_residuals(q[, 4:9], q[, 1:3])
  c <- cor(e)
  setNames(c(c[3, 1], c[4, 1], c[5, 1], c[6, 1], c[3, 2], c[6, 2], c[3, 1] - c[6, 1]),
    pcnames)
}

basecases <- list(pre_full = P0, pre_no2020 = P0[format(P0$Date, "%Y") != "2020", ], pre_2021_22 = P0[P0$Date >=
  as.Date("2021-01-01"), ])

bootrows <- list()
bn <- 0

for (cas in names(basecases)) {
  for (trim in c(FALSE, TRUE)) {
    for (L in c(10, 20)) {
      d0 <- basecases[[cas]]
      v0 <- portfolio_statistics(d0, trim)
      v1 <- portfolio_statistics(P1, trim)
      set.seed(20261014 + L + 100 * trim + match(cas, names(basecases)))
      boot <- replicate(1000, portfolio_statistics(P1[moving_block_indices(nrow(P1),
        L), ], trim) - portfolio_statistics(d0[moving_block_indices(nrow(d0),
        L), ], trim))
      ci <- t(apply(boot, 1, quantile, c(0.025, 0.975)))
      bn <- bn + 1
      bootrows[[bn]] <- data.frame(baseline = cas, winsor1pct = trim, block_days = L,
        metric = pcnames, pre = v0, post = v1, delta = v1 - v0, CI_lower = ci[,
         1], CI_upper = ci[, 2], n_pre = nrow(d0), n_post = nrow(P1), n_boot = 1000)
    }
  }
}

bt <- do.call(rbind, bootrows)

write.csv(bt, file.path(output_dir, "tables", "portfolio_bootstrap_robustness.csv"), row.names = FALSE,
  na = "")

main <- bt[bt$baseline == "pre_full" & !bt$winsor1pct & bt$block_days == 10, ]

write.csv(main, file.path(output_dir, "tables", "primary_partial_correlations.csv"), row.names = FALSE,
  na = "")

contrast <- function(d) {
  e <- factor_residuals(as.matrix(d[, c("chip", "hyper", "merchant", "VST", "NRG", "utility")]),
    d[, c("mkt", "oil", "dy10")])
  c <- cor(e)
  c(merchant_chip_minus_hyper = c[3, 1] - c[3, 2], merchant_chip_minus_utility = c[3,
    1] - c[6, 1], VSTchip_minus_NRGchip = c[4, 1] - c[5, 1])
}

set.seed(5042026)

bb <- replicate(1200, contrast(P1[moving_block_indices(nrow(P1), 10), ]) - contrast(P0[moving_block_indices(nrow(P0),
  10), ]))

ci <- t(apply(bb, 1, quantile, c(0.025, 0.975)))
v0 <- contrast(P0)
v1 <- contrast(P1)

write.csv(data.frame(metric = names(v0), pre = v0, post = v1, change = v1 - v0, CI_low = ci[,
  1], CI_high = ci[, 2]), file.path(output_dir, "tables", "cross_sector_contrast_bootstrap.csv"),
  row.names = FALSE, na = "")

cr <- list()
rr <- list()
ix <- 0
jx <- 0

for (nam in c("pre", "post", "pre_2021_22", "pre_ex2020")) {
  d <- switch(nam, pre = P0, post = P1, pre_2021_22 = basecases$pre_2021_22, pre_ex2020 = basecases$pre_no2020)
  for (ct in list("mkt", c("mkt", "oil"), c("mkt", "oil", "dy10"))) {
    e <- factor_residuals(d[, c("VST", "NRG", "merchant", "utility", "chip", "hyper")],
      d[, ct, drop = FALSE])
    co <- cor(e)
    for (y in 1:4) {
      for (x in 5:6) {
        ix <- ix + 1
        cr[[ix]] <- data.frame(period = nam, controls = paste(ct, collapse = "+"),
         subject = colnames(e)[y], target = colnames(e)[x], partial_corr = co[y,
          x], n_days = nrow(d))
      }
    }
  }
}

write.csv(do.call(rbind, cr), file.path(output_dir, "tables", "controls_partial_correlations.csv"),
  row.names = FALSE, na = "")

panel <- rbind(transform(P0, post = 0), transform(P1, post = 1))

for (resp in c("VST", "NRG", "merchant", "utility")) {
  for (ct in list(c("chip", "mkt"), c("chip", "mkt", "oil"), c("chip", "mkt", "oil",
    "dy10"), c("chip", "hyper", "mkt"), c("chip", "hyper", "mkt", "oil", "dy10"))) {
    f <- lm(as.formula(paste(resp, "~ post * (", paste(ct, collapse = "+"), ")")),
      data = panel)
    ht <- hac_normal_inference(f)
    for (v in intersect(ct, c("chip", "hyper"))) {
      term <- paste0("post:", v)
      row <- ht[ht$term == term, ]
      jx <- jx + 1
      rr[[jx]] <- data.frame(asset = resp, controls = paste(ct, collapse = "+"),
        factor = v, pre_beta = coef(f)[v], post_beta = coef(f)[v] + row$beta,
        delta = row$beta, hac_se = row$se, p = row$p, CI_low = row$low, CI_high = row$high,
        N = nobs(f))
    }
  }
}

reg <- do.call(rbind, rr)
reg$BH_q <- p.adjust(reg$p, "BH")

write.csv(reg, file.path(output_dir, "tables", "HAC_all_interactions.csv"), row.names = FALSE,
  na = "")

write.csv(reg[grepl("hyper", reg$controls), ], file.path(output_dir, "tables", "HAC_chip_vs_hyperscaler_joint.csv"),
  row.names = FALSE, na = "")

write.csv(reg[!grepl("hyper", reg$controls), ], file.path(output_dir, "tables", "HAC_structural_change_exploratory.csv"),
  row.names = FALSE, na = "")

sets <- list(all13 = chip_tickers, big4 = c("NVDA", "AMD", "AVGO", "TSM"), ex_NVDA = setdiff(chip_tickers,
  "NVDA"), ex_SMCI = setdiff(chip_tickers, "SMCI"), classic_semis = c("NVDA", "AMD",
  "AVGO", "INTC", "QCOM", "TXN", "MU", "AMAT", "LRCX", "ADI", "KLAC", "MRVL", "TSM"))

periods <- list(pre_2018_22 = P0, pre_2021_22 = basecases$pre_2021_22, post_2023_26 = P1,
  post_2023_24 = P1[P1$Date < as.Date("2025-01-01"), ], post_2025_26 = P1[P1$Date >=
    as.Date("2025-01-01"), ])

fr <- list()
idx <- 0

for (nam in names(periods)) {
  for (bas in names(sets)) {
    for (ct in list("mkt", c("mkt", "oil", "dy10"))) {
      d <- periods[[nam]]
      yy <- cbind(chip = rowMeans(d[, sets[[bas]]]), hyper = d$hyper, VST = d$VST,
        NRG = d$NRG, utility = d$utility)
      e <- factor_residuals(yy, d[, ct, drop = FALSE])
      co <- cor(e)
      for (j in 3:5) {
        idx <- idx + 1
        fr[[idx]] <- data.frame(period = nam, chip_basket = bas, controls = paste(ct,
         collapse = "+"), firm = colnames(e)[j], chip_r = co[j, 1], hyper_r = co[j,
         2], n = nrow(d))
      }
    }
  }
}

write.csv(do.call(rbind, fr), file.path(output_dir, "tables", "firm_and_basket_robustness.csv"),
  row.names = FALSE, na = "")

group_dependence_statistics <- function(y) {
  co <- cor(y)
  out <- c()
  for (a in names(portfolio_groups)) {
    for (b in names(portfolio_groups)) {
      i <- which(colnames(y) %in% portfolio_groups[[a]])
      j <- which(colnames(y) %in% portfolio_groups[[b]])
      if (a != b && length(i) && length(j)) {
        q <- co[i, j, drop = FALSE]
        key <- paste(a, b, sep = "__")
        out[paste0(key, "__r")] <- mean(q)
        out[paste0(key, "__abs_r")] <- mean(abs(q))
        out[paste0(key, "__r2")] <- mean(q^2)
      }
    }
  }
  theta <- co^2/rowSums(co^2)
  for (a in names(portfolio_groups)) {
    for (b in names(portfolio_groups)) {
      i <- which(colnames(y) %in% portfolio_groups[[a]])
      j <- which(colnames(y) %in% portfolio_groups[[b]])
      if (a != b && length(i) && length(j)) {
        out[paste0(a, "_to_", b, "_GFEVDpct")] <- 100 * mean(rowSums(theta[j,
         i, drop = FALSE]))
      }
    }
  }
  out
}

adjust_dependence_panel <- function(y, m, method) {
  if (method == "raw") {
    return(scale(y, scale = FALSE))
  }
  if (method == "capm") {
    return(factor_residuals(y, m))
  }
  z <- scale(y)
  s <- svd(z, nu = 3, nv = 3)
  z - s$u %*% diag(s$d[1:3]) %*% t(s$v)
}

cases <- list(pre_full = pre_returns, pre_ex2020 = pre_returns[format(as.Date(rownames(pre_returns)),
  "%Y") != "2020", ], pre_2018_19 = pre_returns[rownames(pre_returns) < "2020-01-01",
  ], pre_2021_22 = pre_returns[rownames(pre_returns) >= "2021-01-01", ], post_all = post_returns,
  post_to2025 = post_returns[rownames(post_returns) < "2026-01-01", ])

sr <- list()
idx <- 0
arr <- list()

for (nam in names(cases)) {
  d <- merge(data.frame(Date = as.Date(rownames(cases[[nam]])), cases[[nam]]), market_controls$sp,
    by = "Date")
  d <- d[complete.cases(d), ]
  y <- as.matrix(d[, common])
  arr[[nam]] <- list(y = y, m = d$mkt, date = d$Date)
  for (method in c("raw", "capm", "pca3")) {
    idx <- idx + 1
    sr[[idx]] <- data.frame(period = nam, method = method, N_days = nrow(y), as.list(group_dependence_statistics(adjust_dependence_panel(y,
      d$mkt, method))), check.names = FALSE)
  }
}

ff <- read.csv(file.path(data_dir, "R_FF3_post_ai_returns.csv"), row.names = 1, check.names = FALSE)[,
  common]
idx <- idx + 1

sr[[idx]] <- data.frame(period = "post_all", method = "R_FF3", N_days = nrow(ff), as.list(group_dependence_statistics(as.matrix(ff))),
  check.names = FALSE)

write.csv(do.call(rbind, sr), file.path(output_dir, "tables", "study_all_specifications.csv"),
  row.names = FALSE, na = "")

for (nam in c("pre_full", "pre_ex2020", "pre_2021_22")) {
  d0 <- arr[[nam]]
  d1 <- arr$post_all
  s0 <- group_dependence_statistics(factor_residuals(d0$y, d0$m))
  s1 <- group_dependence_statistics(factor_residuals(d1$y, d1$m))
  keys <- names(s0)[grepl("__r$|GFEVDpct$", names(s0))]
  set.seed(1728)
  boot <- replicate(1200, {
    i <- moving_block_indices(nrow(d0$y), 10)
    j <- moving_block_indices(nrow(d1$y), 10)
    group_dependence_statistics(factor_residuals(d1$y[j, ], d1$m[j]))[keys] - group_dependence_statistics(factor_residuals(d0$y[i,
      ], d0$m[i]))[keys]
  })
  ci <- t(apply(boot, 1, quantile, c(0.025, 0.975)))
  write.csv(data.frame(metric = keys, pre = s0[keys], post = s1[keys], difference = s1[keys] -
    s0[keys], ci_low = ci[, 1], ci_high = ci[, 2], prob_diff_le0 = rowMeans(boot <=
    0), bootstrap_n = 1200, block = 10, seed = 1728), file.path(output_dir, "tables",
    paste0("bootstrap_", nam, ".csv")), row.names = FALSE, na = "")
}

roll <- list()
idx <- 0

for (nam in c("pre_full", "post_all")) {
  d <- arr[[nam]]
  for (e in seq(150, nrow(d$y), by = 21)) {
    ix <- (e - 149):e
    co <- cor(factor_residuals(d$y[ix, ], d$m[ix]))
    idx <- idx + 1
    roll[[idx]] <- data.frame(end = d$date[e], period = nam, merchant__chip__r = mean(co[c("VST",
      "NRG"), chip_tickers]), utilities__chip__r = mean(co[utility_tickers, chip_tickers]))
  }
}

rolling <- do.call(rbind, roll)

write.csv(rolling, file.path(output_dir, "tables", "rolling_150day_sector_exposures.csv"),
  row.names = FALSE, na = "")

# The rolling figure is drawn once below from its exported numerical source.

q <- main[main$metric %in% pcnames[1:4], ]

pdf(file.path(output_dir, "figures", "firm_exposure_changes.pdf"), width = 8, height = 5)
par(mar = c(5, 5, 3, 1), las = 1, cex = 0.85)

bp <- barplot(q$delta, names.arg = q$metric, col = "#457393", ylim = range(c(0, q$CI_lower,
  q$CI_upper)), ylab = "Post minus pre partial correlation")

arrows(bp, q$CI_lower, bp, q$CI_upper, angle = 90, code = 3, length = 0.05)

abline(h = 0)

dev.off()

ffr <- do.call(rbind, lapply(c("VST", "NRG"), function(s) {
  do.call(rbind, lapply(c("chip", "hyperscaler", "Technology", "dc_core"), function(g) {
    members <- intersect(portfolio_groups[[g]], colnames(ff))
    data.frame(asset = s, basket = g, post_FF3_index_correlation = cor(ff[, s], rowMeans(ff[,
      members])), post_FF3_mean_pairwise_corr = mean(cor(ff[, s], ff[, members])))
  }))
}))

write.csv(ffr, file.path(output_dir, "tables", "post_FF3_check.csv"), row.names = FALSE,
  na = "")

lagrows <- list()
idx <- 0

for (nam in c("pre", "post")) {
  d <- if (nam == "pre") {
    P0
  }
  else {
    P1
  }
  e <- factor_residuals(d[, c("merchant", "chip", "utility")], d$mkt)
  l <- data.frame(e[-1, ], lag_merchant = e[-nrow(e), 1], lag_chip = e[-nrow(e), 2],
    lag_utility = e[-nrow(e), 3])
  for (resp in colnames(e)) {
    f <- lm(as.formula(paste(resp, "~ lag_merchant+lag_chip+lag_utility")), data = l)
    ht <- hac_normal_inference(f)
    for (pred in setdiff(colnames(e), resp)) {
      z <- ht[ht$term == paste0("lag_", pred), ]
      idx <- idx + 1
      lagrows[[idx]] <- data.frame(period = nam, target_today = resp, pred_yesterday = pred,
        beta = z$beta, hac_p = z$p, CI_low = z$low, CI_high = z$high, R2 = summary(f)$r.squared,
        N = nobs(f))
    }
  }
}

write.csv(do.call(rbind, lagrows), file.path(output_dir, "tables", "next_day_directional_tests.csv"),
  row.names = FALSE, na = "")

saveRDS(list(P0 = P0, P1 = P1, arr = arr), file.path(output_dir, "models/exposure_inputs.rds"))

# Rolling exposure figure
# 150 observations, 21-day steps; S&P-only residual mean pairwise correlations.

rolling <- read.csv(file.path(output_dir, "tables/rolling_150day_sector_exposures.csv"))

rolling$end <- as.Date(rolling$end)

pdf(file.path(output_dir, "figures", "rolling_merchant_vs_utilities.pdf"), width = 10,
  height = 5.5)
par(mar = c(5, 5, 3, 1), las = 1, cex = 0.85)

plot(range(rolling$end), range(rolling[, 3:4]), type = "n", xlab = "Window end", ylab = "Mean pairwise residual correlation")

rect(as.Date("2022-12-01"), -2, as.Date("2023-09-17"), 2, col = "#eeeeee", border = NA)

for (nam in c("pre_full", "post_all")) {
  q <- rolling[rolling$period == nam, ]
  lines(q$end, q[, 3], col = "#245a81", lwd = 2)
  lines(q$end, q[, 4], col = "#a84b3a", lwd = 2)
}

abline(h = 0, lty = 3)

legend("bottomright", c("VST-NRG / semiconductors", "Utilities / semiconductors"), col = c("#245a81",
  "#a84b3a"), lty = 1, bty = "n")

dev.off()

