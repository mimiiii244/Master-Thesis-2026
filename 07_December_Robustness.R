setwd("/Users/amyhu/Research/AI_Energy_Thesis_Code")
# 07 December Robustness
# Change the working-directory line above for another computer.
data_dir <- "09_Shared_Data/inputs"
output_dir <- "07_December_Robustness/outputs"
for (folder_name in c("tables", "models", "figures", "data")) {
  dir.create(file.path(output_dir, folder_name), recursive = TRUE, showWarnings = FALSE)
}
source("10_Shared_Functions/Shared_Functions.R")
library(BigVAR)

file.copy(file.path(data_dir, "universe.csv"), file.path(output_dir, "tables", "universe.csv"),
  overwrite = TRUE)

# Reconstruct early returns
# Keep the original observed-price construction and the frozen September panel.

options(stringsAsFactors = FALSE)

recovery_tickers <- read.csv(file.path(output_dir, "tables/universe.csv"))$ticker
rawdir <- data_dir
getprice <- function(t) {
  d <- read.csv(file.path(rawdir, paste0("raw_daily_", t, ".csv")), skip = 3, header = FALSE,
    check.names = FALSE)
  data.frame(date = d[, 1], price = as.numeric(d[, 2]))
}

P <- lapply(recovery_tickers, getprice)
names(P) <- recovery_tickers
P <- Map(function(x, n) {
  names(x)[2] <- n
  x
}, P, recovery_tickers)

px <- Reduce(function(a, b) {
  merge(a, b, by = "date", all = TRUE)
}, P)

px <- px[order(px$date), ]
prices <- as.matrix(px[, recovery_tickers])
rownames(prices) <- px$date
r <- log(prices[-1, , drop = FALSE]/prices[-nrow(prices), , drop = FALSE])
rownames(r) <- rownames(prices)[-1]

late <- as.matrix(read.csv(file.path(data_dir, "R_ai_energy_lasso_100_daily_returns.csv"),
  row.names = 1, check.names = FALSE))[, recovery_tickers]

early <- r[rownames(r) >= "2022-12-01" & rownames(r) < "2023-09-15", recovery_tickers,
  drop = FALSE]

stopifnot(nrow(early) >= 165, !anyNA(early), all(is.finite(early)), min(rownames(early)) ==
  "2022-12-01")

ext <- rbind(early, late)

stopifnot(!anyDuplicated(rownames(ext)), all(diff(as.Date(rownames(ext))) > 0))

write.csv(ext, file.path(output_dir, "data", "extended_post90_returns.csv"))

write.csv(data.frame(period = c("early", "post_sep", "post_dec"), n = c(nrow(early), nrow(late),
  nrow(ext)), first = c(min(rownames(early)), min(rownames(late)), min(rownames(ext))),
  last = c(max(rownames(early)), max(rownames(late)), max(rownames(ext)))), file.path(output_dir,
  "tables/recovered_sample_audit.csv"), row.names = FALSE)

write.csv(data.frame(ticker = recovery_tickers, missing_early_returns = colSums(is.na(early))),
  file.path(output_dir, "tables/recovered_early_missingness.csv"), row.names = FALSE)

write.csv(prices, file.path(output_dir, "data", "recovered_frozen_adjusted_prices90.csv"))

# Re-estimate pre and September benchmarks
# The December exercise remains a separate sensitivity.

pre <- as.matrix(read.csv(file.path(data_dir, "R_pre_ai_returns.csv"), row.names = 1,
  check.names = FALSE))
post <- as.matrix(read.csv(file.path(data_dir, "R_ai_energy_lasso_100_daily_returns.csv"),
  row.names = 1, check.names = FALSE))
ff <- read.csv(file.path(data_dir, "FF3_daily_R.csv"))
matched_tickers <- unlist(sector_groups, use.names = FALSE)

matched_tickers <- matched_tickers[matched_tickers %in% intersect(colnames(pre), colnames(post))]

stopifnot(length(matched_tickers) == 90, setequal(matched_tickers, read.csv(file.path(output_dir,
  "tables/universe.csv"))$ticker))

period_fits <- list()
rows <- list()

for (tag in c("pre_main", "post_sep")) {
  r <- if (tag == "pre_main") {
    pre[, matched_tickers, drop = FALSE]
  }
  else {
    post[, matched_tickers, drop = FALSE]
  }
  f <- ff[match(rownames(r), ff$date), ]
  stopifnot(!anyNA(r), !anyNA(f))
  period_fits[[tag]] <- fit_cv_network(adjust_ff_factors(r, f), paste0("baseline_",
    tag), TRUE)
  m <- period_fits[[tag]]
  s <- matched_sector_allocations(m$D)
  z <- matched_sector_allocations(m$Z)
  alloc <- function(x, g) {
    x$allocation[x$receiver == g & x$sender == "Technology"]
  }
  rows[[tag]] <- data.frame(period = tag, n = nrow(r), k = ncol(r), first = min(rownames(r)),
    last = max(rownames(r)), lambda = m$lambda, nonzero = sum(abs(m$Phi) > 1e-12),
    TCI = matched_tci(m$D), generation = alloc(s, "Generation"), utilities = alloc(s,
      "Utilities"), cov_TCI = matched_tci(m$Z), cov_generation = alloc(z, "Generation"),
    cov_utilities = alloc(z, "Utilities"))
}

a <- do.call(rbind, rows)

write.csv(a, file.path(output_dir, "tables", "executed_benchmark_summary.csv"), row.names = FALSE)

saveRDS(period_fits, file.path(output_dir, "models/executed_benchmarks.rds"))

# December-start estimation
# 851 daily observations; period-specific FF3 and fresh 19-point CV.

f <- readRDS(file.path(output_dir, "models/executed_benchmarks.rds"))
r <- as.matrix(read.csv(file.path(file.path(output_dir, "data"), "extended_post90_returns.csv"),
  row.names = 1, check.names = FALSE))
matched_tickers <- colnames(f$pre_main$y)
r <- r[, matched_tickers]
ff <- read.csv(file.path(data_dir, "FF3_daily_R.csv"))
ff <- ff[match(rownames(r), ff$date), ]

stopifnot(!anyNA(ff), nrow(r) == 851, ncol(r) == 90)

f$post_dec <- fit_cv_network(adjust_ff_factors(r, ff), "post_dec_recovered", TRUE)

saveRDS(f, file.path(output_dir, "models/completed_early_start_models.rds"))

saveRDS(list(r = r, f = ff), file.path(output_dir, "models/data_post_dec.rds"))

# December network exports
# all matrices, sector allocations, sparsity and CV boundaries.

period_fits <- readRDS(file.path(output_dir, "models/completed_early_start_models.rds"))

rr <- list(pre_main = as.matrix(read.csv(file.path(data_dir, "R_pre_ai_returns.csv"),
  row.names = 1, check.names = FALSE)), post_sep = as.matrix(read.csv(file.path(data_dir,
  "R_ai_energy_lasso_100_daily_returns.csv"), row.names = 1, check.names = FALSE)),
  post_dec = as.matrix(read.csv(file.path(file.path(output_dir, "data"), "extended_post90_returns.csv"),
    row.names = 1, check.names = FALSE)))

a <- do.call(rbind, lapply(names(period_fits), function(tag) {
  m <- period_fits[[tag]]
  r <- rr[[tag]]
  data.frame(period = tag, n = nrow(r), k = ncol(m$Phi), first = min(rownames(r)), last = max(rownames(r)),
    lambda = m$lambda, nonzero = sum(abs(m$Phi) > 1e-12), sparsity = mean(abs(m$Phi) <=
      1e-12), TCI = matched_tci(m$D), cov_TCI = matched_tci(m$Z), max_cell_diff_pp = 100 *
      max(abs(m$D - m$Z)))
}))

write.csv(a, file.path(output_dir, "tables", "completed_model_summary.csv"), row.names = FALSE)

s <- do.call(rbind, lapply(names(period_fits), function(tag) {
  do.call(rbind, lapply(c("D", "Z"), function(type) {
    x <- matched_sector_allocations(period_fits[[tag]][[type]])
    x$period <- tag
    x$type <- type
    x
  }))
}))

write.csv(s, file.path(output_dir, "tables", "completed_sector_allocations.csv"), row.names = FALSE)

extract_point_metrics <- function(tag, type) {
  x <- s[s$period == tag & s$type == type & s$sender == "Technology", ]
  g <- x$allocation[x$receiver == "Generation"]
  u <- x$allocation[x$receiver == "Utilities"]
  c(generation = g, utilities = u, contrast = g - u)
}

z <- do.call(rbind, lapply(c("D", "Z"), function(type) {
  do.call(rbind, lapply(c("post_sep", "post_dec"), function(tag) {
    p <- extract_point_metrics("pre_main", type)
    v <- extract_point_metrics(tag, type)
    data.frame(post = tag, type = type, generation = v[1], utilities = v[2], contrast = v[3],
      pre_generation = p[1], pre_utilities = p[2], generation_change = v[1] - p[1],
      utility_change = v[2] - p[2], difference_in_changes = v[3] - p[3])
  }))
}))

write.csv(z, file.path(output_dir, "tables", "completed_main_comparison.csv"), row.names = FALSE)

write.csv(data.frame(comparison = "December minus September, common end 2026-04-24", generation = z$generation[2] -
  z$generation[1], utilities = z$utilities[2] - z$utilities[1], difference_in_changes = z$difference_in_changes[2] -
  z$difference_in_changes[1]), file.path(output_dir, "tables", "completed_common_end_comparison.csv"),
  row.names = FALSE)

for (tag in names(period_fits)) {
  for (type in c("D", "Z")) {
    write.csv(period_fits[[tag]][[type]], file.path(output_dir, "tables", paste0("completed_matrix_",
      tag, "_", type, ".csv")), row.names = TRUE)
  }
}

grid <- do.call(rbind, lapply(names(period_fits), function(tag) {
  g <- period_fits[[tag]]$fit@LambdaGrid
  lam <- period_fits[[tag]]$lambda
  data.frame(period = tag, lambda = lam, grid_count = length(g), grid_min = min(g),
    grid_max = max(g), selected_at_lower = abs(lam - min(g)) < 1e-12, selected_at_upper = abs(lam -
      max(g)) < 1e-12)
}))

write.csv(grid, file.path(output_dir, "tables", "completed_CV_grid_diagnostics.csv"),
  row.names = FALSE)

# Early-start uncertainty
# 399 draws, three periods, blocks 10/20: 2,394 refits. Preserve shared pre draws.
# No independent-period interval for the overlapping December–September difference.

period_models <- readRDS(file.path(output_dir, "models/completed_early_start_models.rds"))

B <- 399L

stopifnot(B >= 399)

matched_tickers <- colnames(period_models$pre_main$y)
ff <- read.csv(file.path(data_dir, "FF3_daily_R.csv"))

rr <- list(pre_main = as.matrix(read.csv(file.path(data_dir, "R_pre_ai_returns.csv"),
  row.names = 1, check.names = FALSE))[, matched_tickers], post_sep = as.matrix(read.csv(file.path(data_dir,
  "R_ai_energy_lasso_100_daily_returns.csv"), row.names = 1, check.names = FALSE))[,
  matched_tickers], post_dec = as.matrix(read.csv(file.path(file.path(output_dir, "data"),
  "extended_post90_returns.csv"), row.names = 1, check.names = FALSE))[, matched_tickers])

for (tag in names(rr)) {
  r <- rr[[tag]]
  f <- ff[match(rownames(r), ff$date), ]
  stopifnot(!anyNA(f))
  saveRDS(list(r = r, f = f), paste0(file.path(output_dir, "models/data_"), tag, ".rds"))
}

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
  write.csv(q, file.path(output_dir, "tables", paste0("early_bootstrap_", tag, "_L",
    L, ".csv")), row.names = FALSE)
  q
}

out <- parallel::mclapply(seq_len(nrow(jobs)), run_bootstrap_job, mc.cores = if (.Platform$OS.type ==
  "windows") 1L else 2, mc.preschedule = FALSE)

stopifnot(!any(vapply(out, inherits, logical(1), "try-error")))

q <- do.call(rbind, out)

saveRDS(q, file.path(output_dir, "models/early_bootstrap_draws.rds"))

interval <- function(x) {
  unname(quantile(x, c(0.025, 0.975)))
}

extract_bootstrap_metrics <- function(df, tag, L, type) {
  s <- df[df$period == tag & df$block == L & df$type == type & df$sender == "Technology",
    ]
  g <- s[s$receiver == "Generation", ]
  u <- s[s$receiver == "Utilities", ]
  g <- g[order(g$draw), ]
  u <- u[order(u$draw), ]
  stopifnot(nrow(g) == B, nrow(u) == B, identical(g$draw, u$draw))
  cbind(generation = g$allocation, utilities = u$allocation, contrast = g$allocation -
    u$allocation, TCI = g$TCI)
}

point <- function(tag, type) {
  s <- matched_sector_allocations(period_models[[tag]][[type]])
  g <- s$allocation[s$receiver == "Generation" & s$sender == "Technology"]
  u <- s$allocation[s$receiver == "Utilities" & s$sender == "Technology"]
  c(generation = g, utilities = u, contrast = g - u, TCI = matched_tci(period_models[[tag]][[type]]))
}

ci <- list()
j <- 0

for (L in c(10, 20)) {
  for (type in c("D", "Z")) {
    for (tag in c("post_sep", "post_dec")) {
      a <- extract_bootstrap_metrics(q, "pre_main", L, type)
      b <- extract_bootstrap_metrics(q, tag, L, type)
      delta <- b - a
      p <- point(tag, type) - point("pre_main", type)
      for (v in colnames(delta)) {
        bounds <- interval(delta[, v])
        j <- j + 1
        ci[[j]] <- data.frame(post = tag, block = L, type = type, metric = v,
         point = p[v], bootstrap_mean = mean(delta[, v]), low = bounds[1], high = bounds[2],
         B = B)
      }
    }
  }
}

write.csv(do.call(rbind, ci), file.path(output_dir, "tables", "early_start_bootstrap_change_intervals.csv"),
  row.names = FALSE)

write.csv(aggregate(active ~ period + type + block, q, function(x) {
  c(mean = mean(x), min = min(x), max = max(x))
}), file.path(output_dir, "tables", "early_bootstrap_active_summary.csv"), row.names = FALSE)

