setwd("/Users/amyhu/Research/AI_Energy_Thesis_Code")
# 03 Pre Post Connectedness
# Change the working-directory line above for another computer.
data_dir <- "09_Shared_Data/inputs"
output_dir <- "03_Pre_Post_Connectedness/outputs"
for (folder_name in c("tables", "models", "figures", "data")) {
  dir.create(file.path(output_dir, folder_name), recursive = TRUE, showWarnings = FALSE)
}
source("10_Shared_Functions/Shared_Functions.R")
library(BigVAR)

# Matched periods and fitted networks
# Use the same ordered 90 tickers; period-specific FF3 and fresh CV. Keep the
# archived 94-stock input diagnostic separate.

pre_returns <- as.matrix(read.csv(file.path(data_dir, "R_pre_ai_returns.csv"), row.names = 1,
  check.names = FALSE))
post_returns <- as.matrix(read.csv(file.path(data_dir, "R_ai_energy_lasso_100_daily_returns.csv"),
  row.names = 1, check.names = FALSE))
factor_data <- read.csv(file.path(data_dir, "FF3_daily_R.csv"))
matched_tickers <- unlist(sector_groups, use.names = FALSE)

matched_tickers <- matched_tickers[matched_tickers %in% intersect(colnames(pre_returns),
  colnames(post_returns))]

stopifnot(length(matched_tickers) == 90, setequal(setdiff(unlist(sector_groups), matched_tickers),
  c("CEG", "TLN", "ARM", "PLTR")))

write.csv(data.frame(ticker = matched_tickers, sector = rep(names(sector_groups), lengths(sector_groups))[match(matched_tickers,
  unlist(sector_groups))]), file.path(output_dir, "tables", "universe.csv"), row.names = FALSE)

align_period_factors <- function(r) {
  f <- factor_data[match(rownames(r), factor_data$date), ]
  stopifnot(!anyNA(f), all(is.finite(r)))
  f
}

periods <- list(pre_main = pre_returns[, matched_tickers], post = post_returns[, matched_tickers],
  pre_2018_2019 = pre_returns[substr(rownames(pre_returns), 1, 4) %in% c("2018", "2019"),
    matched_tickers], pre_2021_2022 = pre_returns[substr(rownames(pre_returns), 1,
    4) %in% c("2021", "2022"), matched_tickers])

manifest <- list()
period_models <- list()

for (tag in names(periods)) {
  r <- periods[[tag]]
  f <- align_period_factors(r)
  y <- adjust_ff_factors(r, f)
  period_models[[tag]] <- fit_cv_network(y, tag)
  manifest[[tag]] <- data.frame(period = tag, n = nrow(r), k = ncol(r), start = min(rownames(r)),
    end = max(rownames(r)))
  saveRDS(list(r = r, f = f), paste0(file.path(output_dir, "models/data_"), tag, ".rds"))
}

write.csv(do.call(rbind, manifest), file.path(output_dir, "tables", "samples.csv"), row.names = FALSE)

saveRDS(period_models, file.path(output_dir, "models/main_models.rds"))

# Covariance-only comparisons
# Squared residual correlations, row-normalised; a diagnostic rather than a new causal estimator.

p <- as.matrix(read.csv(file.path(data_dir, "R_pre_ai_returns.csv"), row.names = 1, check.names = FALSE))
r <- as.matrix(read.csv(file.path(data_dir, "R_ai_energy_lasso_100_daily_returns.csv"),
  row.names = 1, check.names = FALSE))
matched_tickers <- read.csv(file.path(output_dir, "tables/universe.csv"))$ticker
factor_data <- read.csv(file.path(data_dir, "FF3_daily_R.csv"))

sets <- list(pre_main = p, pre_2018_2019 = p[substr(rownames(p), 1, 4) %in% c("2018",
  "2019"), ], pre_2021_2022 = p[substr(rownames(p), 1, 4) %in% c("2021", "2022"), ],
  post = r)

q <- lapply(names(sets), function(tag) {
  a <- sets[[tag]][, matched_tickers]
  y <- adjust_ff_factors(a, factor_data[match(rownames(a), factor_data$date), ])
  s <- matched_sector_allocations(covariance_only_fevd(cov(y)))
  s$period <- tag
  s$TCI <- matched_tci(covariance_only_fevd(cov(y)))
  s
})

write.csv(do.call(rbind, q), file.path(output_dir, "tables", "direct_return_covariance.csv"),
  row.names = FALSE)

# Sector, directional and residual diagnostics
# Keep receiver/sender orientation and all original aggregation conventions.

period_models <- readRDS(file.path(output_dir, "models/main_models.rds"))
rows <- list()
ss <- list()
bench <- list()
fd <- list()

for (tag in names(period_models)) {
  o <- period_models[[tag]]
  f <- o$fit
  y <- o$y
  k <- ncol(y)
  rr <- y[-1, ] - sweep(y[-nrow(y), ] %*% t(o$Phi), 2, o$b[, 1], "+")
  eig <- eigen(o$Sigma, symmetric = TRUE, only.values = TRUE)$values
  rows[[tag]] <- data.frame(period = tag, n = nrow(y), k = k, T1 = f@T1, T2 = f@T2,
    window = f@window.size, lambda = o$lambda, initial_boundary = as.numeric(o$initial@OptimalLambda) ==
      max(o$initial@LambdaGrid), expanded_boundary = o$lambda == max(f@LambdaGrid),
    nonzero = sum(abs(o$Phi) > 1e-12), spectral_radius = max(Mod(eigen(o$Phi, only.values = TRUE)$values)),
    Sigma_min_eigen = min(eig), Sigma_condition = max(eig)/min(eig), max_rowsum_error = max(abs(rowSums(o$D) -
      1)), TCI = matched_tci(o$D), zero_TCI = matched_tci(o$Z), max_dynamic_gap_pp = max(abs(o$D -
      o$Z)) * 100, CV_min = min(colMeans(f@InSampMSFE)), OOS_sum_MSE = mean(f@OOSMSFE),
    OOS_per_asset_MSE = mean(f@OOSMSFE)/k, mean_benchmark_sum_MSE = mean(f@MeanMSFE),
    relative_MSE_vs_mean = mean(f@OOSMSFE)/mean(f@MeanMSFE))
  for (type in c("D", "Z")) {
    d <- o[[type]]
    write.csv(d, file.path(output_dir, "tables", paste0("GFEVD_", tag, "_", type,
      ".csv")), row.names = TRUE)
    s <- matched_sector_allocations(d)
    s$period <- tag
    s$type <- type
    ss[[paste(tag, type)]] <- s
    n <- network_positions(d)
    n$period <- tag
    n$type <- type
    fd[[paste(tag, type)]] <- n
  }
  write.csv(o$Phi, file.path(output_dir, "tables", paste0("Phi_", tag, ".csv")), row.names = TRUE)
  write.csv(o$Sigma, file.path(output_dir, "tables", paste0("Sigma_", tag, ".csv")),
    row.names = TRUE)
  C <- cor(y)
  P <- solve(cov(y))
  PC <- -cov2cor(P)
  diag(PC) <- 1
  for (a in names(sector_groups)) {
    for (b in names(sector_groups)) {
      i <- which(colnames(y) %in% sector_groups[[a]])
      j <- which(colnames(y) %in% sector_groups[[b]])
      z <- C[i, j, drop = FALSE]
      p <- PC[i, j, drop = FALSE]
      if (a == b) {
        z <- z[row(z) != col(z)]
        p <- p[row(p) != col(p)]
      }
      bench[[length(bench) + 1]] <- data.frame(period = tag, receiver = a, sender = b,
        mean_corr = mean(z), mean_squared_corr = mean(z^2), mean_partial_corr = mean(p),
        mean_squared_partial = mean(p^2))
    }
  }
  errs <- list()
  for (t in seq(f@T2, nrow(y) - 1)) {
    tr <- y[seq_len(t), , drop = FALSE]
    actual <- y[t + 1, ]
    errs[[length(errs) + 1]] <- data.frame(origin = t, mean_MSE = mean((actual - colMeans(tr))^2),
      zero_MSE = mean(actual^2))
  }
  write.csv(do.call(rbind, errs), file.path(output_dir, "tables", paste0("forecast_benchmarks_",
    tag, ".csv")), row.names = FALSE)
  ac <- sapply(seq_len(k), function(j) {
    acf(rr[, j], lag.max = 5, plot = FALSE)$acf[2:6]
  })
  write.csv(data.frame(ticker = colnames(y), max_abs_acf5 = apply(abs(ac), 2, max),
    LB5_p = sapply(seq_len(k), function(j) {
      Box.test(rr[, j], lag = 5, type = "Ljung-Box")$p.value
    })), file.path(output_dir, "tables", paste0("innovation_serial_", tag, ".csv")),
    row.names = FALSE)
}

diag <- do.call(rbind, rows)

write.csv(diag, file.path(output_dir, "tables", "diagnostics.csv"), row.names = FALSE)

s <- do.call(rbind, ss)

write.csv(s, file.path(output_dir, "tables", "sector_allocations.csv"), row.names = FALSE)

write.csv(do.call(rbind, fd), file.path(output_dir, "tables", "firm_NET.csv"), row.names = FALSE)

write.csv(do.call(rbind, bench), file.path(output_dir, "tables", "correlation_benchmarks.csv"),
  row.names = FALSE)

diffs <- list()

for (pre_returns in setdiff(names(period_models), "post")) {
  for (type in c("D", "Z")) {
    a <- s[s$period == "post" & s$type == type, ]
    b <- s[s$period == pre_returns & s$type == type, ]
    stopifnot(identical(a$receiver, b$receiver) && identical(a$sender, b$sender))
    a$pre <- pre_returns
    a$pre_allocation <- b$allocation
    a$delta_pp <- a$allocation - b$allocation
    diffs[[length(diffs) + 1]] <- a
  }
}

write.csv(do.call(rbind, diffs), file.path(output_dir, "tables", "sector_changes.csv"),
  row.names = FALSE)

h <- list()

for (tag in names(period_models)) {
  for (H in c(1, 5, 10, 20)) {
    o <- period_models[[tag]]
    d <- matched_connectedness(o$Phi, o$Sigma, H)
    h[[length(h) + 1]] <- data.frame(period = tag, H = H, TCI = matched_tci(d), max_gap_from_H1_pp = 100 *
      max(abs(d - matched_connectedness(o$Phi, o$Sigma, 1))))
  }
}

write.csv(do.call(rbind, h), file.path(output_dir, "tables", "horizon_sensitivity.csv"),
  row.names = FALSE)

firmflows <- list()

for (tag in names(period_models)) {
  for (type in c("D", "Z")) {
    d <- period_models[[tag]][[type]]
    matched_tickers <- colnames(d)
    for (g in names(sector_groups)) {
      j <- which(matched_tickers %in% sector_groups[[g]])
      firmflows[[length(firmflows) + 1]] <- data.frame(period = tag, type = type,
        ticker = matched_tickers, sector = g, receives_from_sector_pct = 100 *
         sapply(seq_along(matched_tickers), function(i) {
          sum(d[i, setdiff(j, i)])
         }), sends_to_sector_mean_pct = 100 * sapply(seq_along(matched_tickers),
         function(i) {
          v <- setdiff(j, i)
          if (length(v)) {
           mean(d[v, i])
          }
          else {
           NA_real_
          }
         }))
    }
  }
}

write.csv(do.call(rbind, firmflows), file.path(output_dir, "tables", "firm_sector_allocations.csv"),
  row.names = FALSE)

plateaus <- lapply(names(period_models), function(tag) {
  o <- period_models[[tag]]
  f <- o$fit
  loss <- as.numeric(colMeans(f@InSampMSFE))
  lambda <- as.numeric(f@LambdaGrid)
  best <- min(loss)
  data.frame(period = tag, selected_lambda = o$lambda, largest_tested = max(lambda),
    weakest_lambda_within_1e_12 = min(lambda[loss <= best + 1e-12]), strongest_lambda_within_1e_12 = max(lambda[loss <=
      best + 1e-12]), n_within_1e_12 = sum(loss <= best + 1e-12), loss_at_strongest = loss[which.max(lambda)],
    minimum_loss = best)
})

write.csv(do.call(rbind, plateaus), file.path(output_dir, "tables", "CV_plateaus.csv"),
  row.names = FALSE)

fc <- list()

for (tag in names(period_models)) {
  y <- period_models[[tag]]$y
  C <- cor(y)
  P <- -cov2cor(solve(cov(y)))
  diag(P) <- 1
  matched_tickers <- colnames(y)
  for (g in names(sector_groups)) {
    j <- which(matched_tickers %in% sector_groups[[g]])
    for (i in seq_along(matched_tickers)) {
      jj <- setdiff(j, i)
      fc[[length(fc) + 1]] <- data.frame(period = tag, ticker = matched_tickers[i],
        other_sector = g, n_other = length(jj), mean_signed_corr = mean(C[i, jj]),
        mean_squared_corr = mean(C[i, jj]^2), mean_partial_corr = mean(P[i, jj]))
    }
  }
}

write.csv(do.call(rbind, fc), file.path(output_dir, "tables", "firm_correlation_benchmarks.csv"),
  row.names = FALSE)

# Moving-block uncertainty
# B=399, blocks 10/20, all four periods. Jointly resample returns/factors, refit FF3
# and BigVAR at the period-selected fixed penalty; preserve job order and seeds.

period_models <- readRDS(file.path(output_dir, "models/main_models.rds"))
B <- 399L
out <- list()

jobs <- expand.grid(tag = names(period_models), L = c(10, 20), stringsAsFactors = FALSE)

run_bootstrap_job <- function(j) {
  tag <- jobs$tag[j]
  L <- jobs$L[j]
  d <- readRDS(paste0(file.path(output_dir, "models/data_"), tag, ".rds"))
  lambda <- period_models[[tag]]$lambda
  set.seed(20261016 + match(tag, names(period_models)) * 100 + L)
  rows <- vector("list", B)
  for (b in seq_len(B)) {
    ix <- moving_block_indices(nrow(d$r), L)
    y <- adjust_ff_factors(d$r[ix, , drop = FALSE], d$f[ix, , drop = FALSE])
    o <- fit_fixed_penalty(y, lambda)
    rows[[b]] <- do.call(rbind, lapply(c("D", "Z"), function(type) {
      s <- matched_sector_allocations(o[[type]])
      s$period <- tag
      s$type <- type
      s$block <- L
      s$draw <- b
      s$TCI <- matched_tci(o[[type]])
      s$active <- sum(abs(o$Phi) > 1e-12)
      s
    }))
  }
  q <- do.call(rbind, rows)
  write.csv(q, file.path(output_dir, "tables", paste0("bootstrap_", tag, "_L", L, ".csv")),
    row.names = FALSE)
  q
}

# Run in the current R session: forked workers can disappear in macOS RStudio.
# Each job sets its original seed, so sequential execution preserves its draws.
out <- lapply(seq_len(nrow(jobs)), function(j) {
  message("Bootstrap: ", jobs$tag[j], ", block ", jobs$L[j])
  result <- run_bootstrap_job(j)
  required <- c("period", "type", "block", "draw", "receiver", "sender",
    "allocation", "TCI", "active")
  stopifnot(is.data.frame(result), all(required %in% names(result)),
    nrow(result) == B * 2L * length(sector_groups)^2,
    all(is.finite(result$allocation)), all(is.finite(result$TCI)),
    all(is.finite(result$active)),
    all(result$period == jobs$tag[j]), all(result$block == jobs$L[j]))
  result
})

q <- do.call(rbind, out)

saveRDS(q, file.path(output_dir, "models/bootstrap_draws.rds"))

ci <- list()

for (pre_returns in setdiff(names(period_models), "post")) {
  for (L in c(10, 20)) {
    for (type in c("D", "Z")) {
      for (a in names(sector_groups)) {
        for (b in names(sector_groups)) {
         p <- q[q$period == "post" & q$block == L & q$type == type & q$receiver ==
          a & q$sender == b, ]
         r <- q[q$period == pre_returns & q$block == L & q$type == type & q$receiver ==
          a & q$sender == b, ]
         stopifnot(nrow(p) == B, nrow(r) == B)
         delta <- p$allocation - r$allocation
         point <- matched_sector_allocations(period_models$post[[type]])$allocation[matched_sector_allocations(period_models$post[[type]])$receiver ==
          a & matched_sector_allocations(period_models$post[[type]])$sender ==
          b] - matched_sector_allocations(period_models[[pre_returns]][[type]])$allocation[matched_sector_allocations(period_models[[pre_returns]][[type]])$receiver ==
          a & matched_sector_allocations(period_models[[pre_returns]][[type]])$sender ==
          b]
         ci[[length(ci) + 1]] <- data.frame(pre = pre_returns, block = L, type = type,
          receiver = a, sender = b, point_delta = point, bootstrap_mean = mean(delta),
          low = quantile(delta, 0.025), high = quantile(delta, 0.975), B = B)
        }
      }
    }
  }
}

write.csv(do.call(rbind, ci), file.path(output_dir, "tables", "bootstrap_change_intervals.csv"),
  row.names = FALSE)

write.csv(aggregate(active ~ period + type + block, q, function(x) {
  c(mean = mean(x), min = min(x), max = max(x))
}), file.path(output_dir, "tables", "bootstrap_active_summary.csv"), row.names = FALSE)

# Fresh-CV composition and factor alternatives
# raw, CAPM, generation-firm exclusion and equal-length pre comparisons.

pre_returns <- as.matrix(read.csv(file.path(data_dir, "R_pre_ai_returns.csv"), row.names = 1,
  check.names = FALSE))
post_returns <- as.matrix(read.csv(file.path(data_dir, "R_ai_energy_lasso_100_daily_returns.csv"),
  row.names = 1, check.names = FALSE))
factor_data <- read.csv(file.path(data_dir, "FF3_daily_R.csv"))
matched_tickers <- read.csv(file.path(output_dir, "tables/universe.csv"))$ticker

periods <- list(pre_main = pre_returns[, matched_tickers], post = post_returns[, matched_tickers])

jobs <- list()

for (tag in names(periods)) {
  for (v in c("raw", "CAPM", "drop_VST", "drop_NRG", "drop_EXC")) {
    jobs[[length(jobs) + 1]] <- list(tag = tag, v = v)
  }
}

jobs[[length(jobs) + 1]] <- list(tag = "pre_equal_n", v = "equal_n_last_pre")
run <- function(job) {
  tag <- job$tag
  v <- job$v
  r <- if (tag == "pre_equal_n") {
    tail(periods$pre_main, nrow(post_returns))
  }
  else {
    periods[[tag]]
  }
  f <- factor_data[match(rownames(r), factor_data$date), ]
  if (grepl("^drop_", v)) {
    r <- r[, setdiff(colnames(r), sub("drop_", "", v)), drop = FALSE]
  }
  mode <- if (v %in% c("raw", "CAPM")) {
    v
  }
  else {
    "FF3"
  }
  o <- fit_cv_network(adjust_ff_factors(r, f, mode), paste0(tag, "_", v))
  s <- matched_sector_allocations(o$D)
  s$period <- tag
  s$variant <- v
  s$TCI <- matched_tci(o$D)
  s$nonzero <- sum(abs(o$Phi) > 1e-12)
  s$start <- min(rownames(r))
  s$end <- max(rownames(r))
  s$n <- nrow(r)
  s
}

res <- lapply(jobs, function(job) {
  message("Fresh CV: ", job$tag, ", ", job$v)
  result <- run(job)
  stopifnot(is.data.frame(result), nrow(result) == length(sector_groups)^2,
    all(is.finite(result$allocation)), all(is.finite(result$TCI)))
  result
})

write.csv(do.call(rbind, res), file.path(output_dir, "tables", "sensitivity_fresh_CV.csv"),
  row.names = FALSE)

# Conditional penalties and asset-count diagnostics
# Lambda multiples .5/.8/1/1.2/2; launch-day exclusion; subset/shift seed 1716.

period_models <- readRDS(file.path(output_dir, "models/main_models.rds"))
sens <- list()
record <- function(o, tag, variant) {
  s <- matched_sector_allocations(o$D)
  s$period <- tag
  s$variant <- variant
  s$TCI <- matched_tci(o$D)
  s$nonzero <- sum(abs(o$Phi) > 1e-12)
  s
}

for (tag in c("pre_main", "post")) {
  for (mult in c(0.5, 0.8, 1, 1.2, 2)) {
    o <- fit_fixed_penalty(period_models[[tag]]$y, period_models[[tag]]$lambda * mult)
    sens[[length(sens) + 1]] <- record(o, tag, paste0("lambda_x", mult))
  }
}

d <- readRDS(file.path(output_dir, "models/data_pre_main.rds"))
i <- rownames(d$r) < "2022-11-30"

o <- fit_fixed_penalty(adjust_ff_factors(d$r[i, ], d$f[i, ]), period_models$pre_main$lambda)

sens[[length(sens) + 1]] <- record(o, "pre_main", "exclude_2022_11_30_fixed_lambda")

write.csv(do.call(rbind, sens), file.path(output_dir, "tables", "sensitivity_sector.csv"),
  row.names = FALSE)

set.seed(1716)

counts <- list()

for (tag in names(period_models)) {
  y <- period_models[[tag]]$y
  for (k in c(10, 20, 40, 60, 90)) {
    for (b in 1:100) {
      j <- sample.int(ncol(y), k)
      z <- covariance_only_fevd(cov(y[, j]))
      counts[[length(counts) + 1]] <- data.frame(period = tag, k = k, draw = b,
        type = "observed_subset", TCI = matched_tci(z))
      yy <- sapply(j, function(j) {
        shift <- sample.int(nrow(y) - 1, 1)
        y[c((shift + 1):nrow(y), 1:shift), j]
      })
      counts[[length(counts) + 1]] <- data.frame(period = tag, k = k, draw = b,
        type = "independent_circular_shifts", TCI = matched_tci(covariance_only_fevd(cov(yy))))
    }
  }
}

write.csv(do.call(rbind, counts), file.path(output_dir, "tables", "asset_count_TCI.csv"),
  row.names = FALSE)

write.csv(aggregate(TCI ~ period + k + type, do.call(rbind, counts), mean), file.path(output_dir,
  "tables", "asset_count_means.csv"), row.names = FALSE)

# Matched-network figures
# Use exported full-precision numerical results; preserve plotting specifications.

pdf <- function(file, width = 7, height = 7, ...) {
  grDevices::cairo_pdf(file, width = width, height = height, family = "DejaVu Sans",
    onefile = TRUE, ...)
}

s <- read.csv(file.path(output_dir, "tables/sector_allocations.csv"))
d <- read.csv(file.path(output_dir, "tables/diagnostics.csv"))
q <- read.csv(file.path(output_dir, "tables/bootstrap_change_intervals.csv"))

pdf(file.path(output_dir, "figures/TCI_dynamic_vs_covariance.pdf"), 9, 5)

par(mar = c(7, 5, 2, 1))

barplot(t(as.matrix(d[, c("TCI", "zero_TCI")])), beside = TRUE, names.arg = c("Main pre",
  "Post", "2018-2019", "2021-Nov 2022"), las = 2, col = c("#226699", "#cc8844"), ylab = "TCI (%)",
  legend.text = c("BigVAR H=10", "Same-Sigma zero lag"))

dev.off()

pdf(file.path(output_dir, "figures/sector_change_heatmaps.pdf"), 9, 8)

for (pre_returns in c("pre_main", "pre_2018_2019", "pre_2021_2022")) {
  a <- s[s$period == "post" & s$type == "D", ]
  b <- s[s$period == pre_returns & s$type == "D", ]
  mat <- matrix(a$allocation - b$allocation, 6, 6, byrow = TRUE)
  lim <- max(abs(mat))
  par(mar = c(9, 9, 4, 2))
  image(1:6, 1:6, t(mat[6:1, ]), col = colorRampPalette(c("#236099", "white", "#b5413d"))(101),
    zlim = c(-lim, lim), axes = FALSE, xlab = "Sender", ylab = "", main = paste("Post minus",
      pre_returns, "(pp); receiver rows"))
  axis(1, 1:6, names(sector_groups), las = 2)
  axis(2, 1:6, rev(names(sector_groups)), las = 2)
  for (i in 1:6) {
    for (j in 1:6) {
      text(j, 7 - i, sprintf("%.2f", mat[i, j]), cex = 0.8)
    }
  }
}

dev.off()

z <- q[q$type == "D" & q$receiver == "Generation" & q$sender == "Technology", ]

pdf(file.path(output_dir, "figures/technology_to_generation_intervals.pdf"), 10, 5)

par(mar = c(5, 11, 4, 2))

labels <- paste(c("Main pre", "Main pre", "2018-2019", "2018-2019", "2021-Nov 2022", "2021-Nov 2022"),
  paste0("(L=", z$block, ")"))

yy <- rev(seq_len(nrow(z)))

plot(z$point_delta, yy, xlim = range(z$low, z$high, 0), ylim = c(0.5, nrow(z) + 0.5),
  yaxt = "n", xlab = "Post minus pre allocation (percentage points)", ylab = "", pch = 19,
  col = "#226699", main = "Technology to Generation: 95% bootstrap intervals")

arrows(z$low, yy, z$high, yy, angle = 90, code = 3, length = 0.05)

axis(2, yy, labels, las = 1)

abline(v = 0, lty = 2)

mtext("399 draws per period; fixed selected penalties", side = 3, line = 0.3, cex = 0.8)

dev.off()

z <- read.csv(file.path(output_dir, "tables/asset_count_means.csv"))

pdf(file.path(output_dir, "figures/asset_count_effect.pdf"), 8, 5)

plot(NA, xlim = c(10, 90), ylim = range(z$TCI), xlab = "Number of assets", ylab = "Covariance-only TCI (%)")

col <- c("#226699", "#cc8844", "#53845c", "#a74563")

for (i in seq_along(unique(z$period))) {
  for (type in unique(z$type)) {
    a <- z[z$period == unique(z$period)[i] & z$type == type, ]
    a <- a[order(a$k), ]
    lines(a$k, a$TCI, col = col[i], lty = if (type == "observed_subset") {
      1
    }
    else {
      2
    })
  }
}

legend("topleft", unique(z$period), col = col, lty = 1, bty = "n")

legend("bottomright", c("Observed", "Independent time shifts"), lty = 1:2, bty = "n")

dev.off()

f <- read.csv(file.path(output_dir, "tables/firm_correlation_benchmarks.csv"))

f <- f[f$other_sector == "Technology" & f$ticker %in% c("VST", "NRG", "EXC") & f$period %in%
  c("pre_main", "post"), ]

m <- sapply(c("pre_main", "post"), function(tag) {
  z <- f[f$period == tag, ]
  z$mean_signed_corr[match(c("VST", "NRG", "EXC"), z$ticker)]
})

pdf(file.path(output_dir, "figures/firm_technology_signed_correlations.pdf"), 8, 5)

barplot(t(m), beside = TRUE, names.arg = c("VST", "NRG", "EXC"), col = c("#226699", "#cc8844"),
  legend.text = c("Main pre", "Post"), ylab = "Mean signed correlation with 29 technology firms",
  main = "Period-specific FF3 residual correlations")

abline(h = 0)

dev.off()

pdf(file.path(output_dir, "figures/CV_loss_paths.pdf"), 10, 8)

par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))

for (tag in d$period) {
  p <- read.csv(paste0(file.path(output_dir, "tables/CV_"), tag, ".csv"))
  plot(p$lambda, 100 * (p$loss/min(p$loss) - 1), log = "x", type = "b", pch = 16, col = "#226699",
    xlab = "Lambda (log scale)", ylab = "Loss above minimum (%)", main = tag)
  points(p$lambda[p$selected], 100 * (p$loss[p$selected]/min(p$loss) - 1), pch = 19,
    col = "#bb4433", cex = 1.5)
}

dev.off()

