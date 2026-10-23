setwd("/Users/amyhu/Research/AI_Energy_Thesis_Code")
# 02 Full Sample Connectedness
# Change the working-directory line above for another computer.
data_dir <- "09_Shared_Data/inputs"
output_dir <- "02_Full_Sample_Connectedness/outputs"
for (folder_name in c("tables", "models", "figures", "data")) {
  dir.create(file.path(output_dir, folder_name), recursive = TRUE, showWarnings = FALSE)
}
source("10_Shared_Functions/Shared_Functions.R")
library(BigVAR)

# Fit 94-stock networks and specification checks
# BigVAR 1.1.5; original grids, lag blocks, intercept, CV timing and seed 42.
# GFEVD sums A0 through A9 for H=10 and normalises receiver rows.

ticker_order <- unlist(sector_groups, use.names = FALSE)

ff3_returns <- as.matrix(read.csv(file.path(data_dir, "R_FF3_ai_energy_lasso_100_daily_returns.csv"),
  row.names = 1, check.names = FALSE))

ff3_returns <- ff3_returns[, ticker_order[ticker_order %in% colnames(ff3_returns)]]
fit_full_network <- function(y, p = 1, struct = "Basic", tag = "baseline") {
  set.seed(42)
  m <- constructModel(y, p = p, struct = struct, gran = c(150, 10), verbose = FALSE,
    IC = tag %in% c("original_inputs", "R_FF3"), VARX = list(), cv = "Rolling", model.control = list(intercept = TRUE))
  modelpath <- file.path(output_dir, "models", paste0("bigvar_", tag, ".rds"))
  f <- cv.BigVAR(m)
  saveRDS(f, modelpath)
  write.csv(data.frame(lambda = as.numeric(f@LambdaGrid), validation_loss = as.numeric(apply(f@InSampMSFE,
    2, mean)), selected = seq_along(f@LambdaGrid) %in% f@index), file.path(output_dir,
    "tables", paste0("CV_path_", tag, ".csv")), row.names = FALSE, na = "")
  write.csv(data.frame(tag = tag, version = as.character(packageVersion("BigVAR")),
    T = nrow(y), k = ncol(y), T1 = f@T1, T2 = f@T2, lambda = as.numeric(f@OptimalLambda),
    index = f@index, grid_n = length(f@LambdaGrid), cv = f@crossval, window = f@window.size,
    h = f@horizon, ONESE = f@ONESE), file.path(output_dir, "tables", paste0("CV_design_",
    tag, ".csv")), row.names = FALSE, na = "")
  b <- f@betaPred
  k <- ncol(y)
  stopifnot(ncol(b) == 1 + k * p)
  phi <- b[, -1, drop = FALSE]
  dimnames(phi) <- list(colnames(y), if (p == 1) {
    colnames(y)
  } else {
    unlist(lapply(seq_len(p), function(j) {
      paste0(colnames(y), "_L", j)
    }))
  })
  x <- lagged_predictors(y, p)
  z <- y[(p + 1):nrow(y), , drop = FALSE]
  r <- sweep(z - x %*% t(phi), 2, b[, 1], "-")
  sig <- crossprod(r)/nrow(r)
  d <- core_connectedness(phi, sig, p)
  dimnames(d) <- list(colnames(y), colnames(y))
  write.csv(phi, file.path(output_dir, "tables", paste0("Phi_", tag, ".csv")), row.names = TRUE,
    na = "")
  write.csv(d, file.path(output_dir, "tables", paste0("GFEVD_", tag, ".csv")), row.names = TRUE,
    na = "")
  write.csv(directional_positions(d, colnames(y)), file.path(output_dir, "tables", paste0("NET_",
    tag, ".csv")), row.names = FALSE, na = "")
  list(fit = f, Phi = phi, Sigma = sig, D = d, y = y, p = p)
}


initial_network <- fit_full_network(ff3_returns, tag = "R_FF3")

saveRDS(initial_network, file.path(output_dir, "models/network_baseline.rds"))

write.csv(core_sector_allocations(initial_network$D, colnames(ff3_returns)), file.path(output_dir,
  "tables", "full94_sector_destination.csv"), row.names = FALSE, na = "")

write.csv(directional_positions(initial_network$D, colnames(ff3_returns)), file.path(output_dir,
  "tables", "full94_net_positions.csv"), row.names = FALSE, na = "")

q <- core_sector_allocations(initial_network$D, colnames(ff3_returns))
q <- q[q$sender == "Utilities", ]

write.csv(q, file.path(output_dir, "tables", "utilities_destinations.csv"), row.names = FALSE,
  na = "")

# The final selected-grid utility destination figure is drawn after penalty selection.

q <- directional_positions(initial_network$D, colnames(ff3_returns))
q <- q[order(q$NET), ]

pdf(file.path(output_dir, "figures", "net_spillovers_initial_grid.pdf"), width = 12, height = 6)
par(mar = c(5, 5, 3, 1), las = 1, cex = 0.85)

par(mar = c(7, 5, 3, 1), las = 2)

barplot(q$NET, names.arg = q$ticker, col = ifelse(q$NET > 0, "#457393", "#ad5948"), ylab = "Net generalized connectedness (percentage points)",
  cex.names = 0.55)

dev.off()

rows <- list()
idx <- 0
add <- function(obj, tag, H = 10) {
  d <- core_connectedness(obj$Phi, obj$Sigma, obj$p, H)
  data.frame(specification = tag, p = obj$p, H = H, nonzero = sum(abs(obj$Phi) > 1e-12),
    TCI = core_tci(d), utility_NET = sum(directional_positions(d, colnames(obj$y))$NET[colnames(obj$y) %in%
      utility_tickers]), max_H1_gap = max(abs(d - core_connectedness(obj$Phi, obj$Sigma,
      obj$p, 1))) * 100)
}

for (H in c(1, 5, 10, 20)) {
  idx <- idx + 1
  rows[[idx]] <- add(initial_network, paste0("baseline_H", H), H)
}

for (p in c(2, 3)) {
  obj <- fit_full_network(ff3_returns, p = p, tag = paste0("p", p))
  idx <- idx + 1
  rows[[idx]] <- add(obj, paste0("p", p))
}

obj <- fit_full_network(ff3_returns, struct = "HLAGC", tag = "HLAGC")
idx <- idx + 1
rows[[idx]] <- add(obj, "HLAGC")

raw <- as.matrix(read.csv(file.path(data_dir, "R_ai_energy_lasso_100_daily_returns.csv"),
  row.names = 1, check.names = FALSE))[rownames(ff3_returns), colnames(ff3_returns)]

mc <- readRDS("09_Shared_Data/frozen_models/market_controls.rds")$sp
mkt <- mc$mkt[match(as.Date(rownames(raw)), mc$Date)]

stopifnot(all(is.finite(mkt)))

for (nam in c("raw", "CAPM")) {
  yy <- if (nam == "raw") {
    raw
  }
  else {
    factor_residuals(raw, mkt)
  }
  dimnames(yy) <- dimnames(raw)
  obj <- fit_full_network(yy, tag = nam)
  idx <- idx + 1
  rows[[idx]] <- add(obj, nam)
}

rob <- do.call(rbind, rows)

write.csv(rob, file.path(output_dir, "tables", "network_robustness.csv"), row.names = FALSE,
  na = "")

pdf(file.path(output_dir, "figures", "network_robustness_initial_grid.pdf"), width = 10,
  height = 5)
par(mar = c(5, 5, 3, 1), las = 1, cex = 0.85)

par(mar = c(7, 5, 3, 1), las = 2)

barplot(rob$TCI, names.arg = rob$specification, col = "#457393", ylab = "TCI (percent)")

dev.off()

post <- read.csv(file.path(data_dir, "R_post_ai_returns.csv"), row.names = 1, check.names = FALSE)

s16 <- c("NVDA", "MSFT", "GOOGL", "META", "AAPL", "AVGO", "AMD", "CEG", "VST", "TLN",
  "NEE", "DUK", "SO", "XOM", "CVX", "COP")

o16 <- fit_ols_var(as.matrix(post[, s16]), 2)
d16 <- core_connectedness(o16$Phi, o16$Sigma, 2)

write.csv(d16, file.path(output_dir, "tables", "stage2_16var_GFEVD.csv"), row.names = TRUE,
  na = "")

write.csv(directional_positions(d16, s16), file.path(output_dir, "tables", "stage2_16var_NET.csv"),
  row.names = FALSE, na = "")

loadseries <- function(f, n, y = FALSE) {
  d <- if (y) {
    read_adjusted_prices(f)
  }
  else {
    read_fred_series(f)
  }
  names(d)[2] <- n
  d[is.finite(d[[2]]), ]
}

bl <- Reduce(function(a, b) {
  merge(a, b, by = "Date")
}, list(loadseries("raw_daily_nasdaq100.csv", "NASDAQ"), loadseries("raw_FRED_broad_usd.csv",
  "USD"), loadseries("raw_dailyish_vix.csv", "VIX"), loadseries("raw_DGS10_daily.csv",
  "DGS10"), loadseries("raw_ai_hardware_smh_daily.csv", "SMH", TRUE), loadseries("raw_ai_hardware_xlu_daily.csv",
  "XLU", TRUE)))

bl <- bl[order(bl$Date), ]
br <- sapply(bl[, -1], diff)

for (n in c("NASDAQ", "USD", "SMH", "XLU")) {
  br[, n] <- diff(log(bl[[n]]))
}

rownames(br) <- as.character(bl$Date[-1])
br <- br[rownames(br) >= "2023-09-15" & rownames(br) <= "2026-04-24", ]
ob <- fit_ols_var(br, 2)
db <- core_connectedness(ob$Phi, ob$Sigma, 2)

write.csv(db, file.path(output_dir, "tables", "stage1_GFEVD.csv"), row.names = TRUE, na = "")

write.csv(directional_positions(db, colnames(br)), file.path(output_dir, "tables", "stage1_NET.csv"),
  row.names = FALSE, na = "")

write.csv(data.frame(model = c("6 financial proxies VAR2", "16 equities VAR2", "94 FF3 BigVAR1"),
  TCI = c(core_tci(db), core_tci(d16), core_tci(initial_network$D))), file.path(output_dir,
  "tables", "model_ladder.csv"), row.names = FALSE, na = "")

# Expanded and refined penalty search
# original candidates, 13-point local refinement and nine upper multiples.

base <- readRDS(file.path(output_dir, "models/network_baseline.rds"))
original <- readRDS(file.path(output_dir, "models/bigvar_R_FF3.rds"))
y <- original@Data
original_penalties <- as.numeric(original@LambdaGrid)
largest_original_penalty <- max(original_penalties)

grid <- sort(unique(c(original_penalties, exp(seq(log(sort(original_penalties, decreasing = TRUE)[2]),
  log(largest_original_penalty), length.out = 13)), largest_original_penalty * c(1.05,
  1.1, 1.25, 1.5, 2, 4, 8, 16, 32))), decreasing = TRUE)

set.seed(42)

model <- constructModel(y, p = 1, struct = "Basic", gran = grid, ownlambdas = TRUE, h = original@horizon,
  cv = original@crossval, T1 = original@T1, T2 = original@T2, window.size = original@window.size,
  ONESE = original@ONESE, verbose = FALSE, IC = FALSE, VARX = list(), model.control = list(intercept = TRUE))

fit <- cv.BigVAR(model)

saveRDS(fit, file.path(output_dir, "models/bigvar_expanded_grid.rds"))

loss <- colMeans(fit@InSampMSFE)

stopifnot(all(original_penalties %in% as.numeric(fit@LambdaGrid)))

path <- data.frame(lambda = as.numeric(fit@LambdaGrid), validation_loss = loss, original_grid = as.numeric(fit@LambdaGrid) %in%
  original_penalties, selected = seq_along(loss) %in% fit@index)

write.csv(path, file.path(output_dir, "tables", "CV_path_expanded_grid.csv"), row.names = FALSE,
  na = "")

write.csv(data.frame(origin = rep(seq_len(nrow(fit@InSampMSFE)), ncol(fit@InSampMSFE)),
  lambda = rep(as.numeric(fit@LambdaGrid), each = nrow(fit@InSampMSFE)), squared_error = as.numeric(fit@InSampMSFE)),
  file.path(output_dir, "tables", "CV_losses_expanded_grid_by_origin.csv"), row.names = FALSE,
  na = "")

path$active_final <- NA_integer_

for (j in seq_len(nrow(path))) {
  fixed <- constructModel(y, p = 1, struct = "Basic", gran = path$lambda[j], ownlambdas = TRUE,
    cv = "None", verbose = FALSE, IC = FALSE, VARX = list(), model.control = list(intercept = TRUE))
  bb <- BigVAR.est(fixed)$B[, , 1]
  path$active_final[j] <- sum(abs(bb[, -1, drop = FALSE]) > 1e-12)
}

write.csv(path, file.path(output_dir, "tables", "CV_path_expanded_grid.csv"), row.names = FALSE,
  na = "")

obj <- function(f) {
  b <- f@betaPred
  phi <- b[, -1, drop = FALSE]
  dimnames(phi) <- list(colnames(y), colnames(y))
  r <- y[-1, ] - sweep(y[-nrow(y), ] %*% t(phi), 2, b[, 1], "+")
  sig <- crossprod(r)/nrow(r)
  d <- core_connectedness(phi, sig)
  dimnames(d) <- list(colnames(y), colnames(y))
  list(fit = f, Phi = phi, Sigma = sig, D = d, y = y, p = 1)
}

base <- obj(original)
new <- obj(fit)

saveRDS(new, file.path(output_dir, "models/network_expanded_grid.rds"))

compare <- function(o, label) {
  d <- o$D
  mm <- directional_positions(d, colnames(y))
  data.frame(specification = label, lambda = o$fit@OptimalLambda, validation_loss = min(colMeans(o$fit@InSampMSFE)),
    evaluation_mse = mean(o$fit@OOSMSFE), active_slopes = sum(abs(o$Phi) > 1e-12),
    TCI = core_tci(d), utility_NET = sum(mm$NET[mm$ticker %in% utility_tickers]),
    max_GFEVD_change_pp = 100 * max(abs(d - base$D)), max_NET_change_pp = max(abs(mm$NET -
      directional_positions(base$D, colnames(y))$NET)))
}

comparison <- rbind(compare(base, "Original ten-point grid"), compare(new, "Expanded and refined grid"))

zero <- new
zero$Phi[, ] <- 0
zero$Sigma <- cov(y[-1, ]) * (nrow(y) - 2)/(nrow(y) - 1)
zero$D <- core_connectedness(zero$Phi, zero$Sigma)
dimnames(zero$D) <- list(colnames(y), colnames(y))

write.csv(zero$D, file.path(output_dir, "tables", "GFEVD_zero_lag.csv"), row.names = TRUE,
  na = "")

write.csv(comparison, file.path(output_dir, "tables", "penalty_boundary_comparison.csv"),
  row.names = FALSE, na = "")

shared <- match(original_penalties, path$lambda)

write.csv(data.frame(lambda = original_penalties, original_loss = colMeans(original@InSampMSFE),
  expanded_loss = path$validation_loss[shared], difference = path$validation_loss[shared] -
    colMeans(original@InSampMSFE)), file.path(output_dir, "tables", "penalty_original_grid_recheck.csv"),
  row.names = FALSE, na = "")

# Selected network and independent zero-lag checks
# Use the expanded fit for main network outputs. Verify intercept-only validation
# separately; evaluation observations do not select the penalty.

extended <- readRDS(file.path(output_dir, "models/network_expanded_grid.rds"))

saveRDS(extended, file.path(output_dir, "models/network_baseline.rds"))

nm <- colnames(extended$y)

write.csv(extended$Phi, file.path(output_dir, "tables", "Phi_expanded_grid.csv"), row.names = TRUE,
  na = "")

write.csv(extended$D, file.path(output_dir, "tables", "GFEVD_expanded_grid.csv"), row.names = TRUE,
  na = "")

write.csv(directional_positions(extended$D, nm), file.path(output_dir, "tables", "NET_expanded_grid.csv"),
  row.names = FALSE, na = "")

write.csv(directional_positions(extended$D, nm), file.path(output_dir, "tables", "full94_net_positions.csv"),
  row.names = FALSE, na = "")

st <- core_sector_allocations(extended$D, nm)

write.csv(st, file.path(output_dir, "tables", "full94_sector_destination.csv"), row.names = FALSE,
  na = "")

q <- st[st$sender == "Utilities", ]

write.csv(q, file.path(output_dir, "tables", "utilities_destinations.csv"), row.names = FALSE,
  na = "")

pdf(file.path(output_dir, "figures", "utility_destinations.pdf"), width = 9, height = 5)
par(mar = c(5, 5, 3, 1), las = 1, cex = 0.85)

q <- q[q$receiver != "Utilities", ]

barplot(q$mean_receiving_pct, names.arg = q$receiver, horiz = TRUE, las = 1, col = "#457393",
  xlab = "Mean percent of receiver variance", cex.names = 0.8)

dev.off()

rob <- read.csv(file.path(output_dir, "tables/network_robustness.csv"))

for (j in which(grepl("^baseline_H", rob$specification))) {
  dd <- core_connectedness(extended$Phi, extended$Sigma, H = rob$H[j])
  mm <- directional_positions(dd, nm)
  rob$nonzero[j] <- sum(abs(extended$Phi) > 1e-12)
  rob$TCI[j] <- core_tci(dd)
  rob$utility_NET[j] <- sum(mm$NET[mm$ticker %in% utility_tickers])
  rob$max_H1_gap[j] <- 100 * max(abs(dd - core_connectedness(extended$Phi, extended$Sigma,
    H = 1)))
}

write.csv(rob, file.path(output_dir, "tables", "network_robustness.csv"), row.names = FALSE,
  na = "")

lad <- read.csv(file.path(output_dir, "tables/model_ladder.csv"))
lad$TCI[3] <- core_tci(extended$D)

write.csv(lad, file.path(output_dir, "tables", "model_ladder.csv"), row.names = FALSE,
  na = "")

ff <- readRDS(file.path(output_dir, "models/bigvar_R_FF3.rds"))
yy <- ff@Data
v <- seq(ff@T1 - 1, ff@T2 - 1)
mean_loss <- vapply(v, function(i) {
  sum((yy[i + 1, ] - colMeans(yy[2:i, , drop = FALSE]))^2)
}, numeric(1))

path <- read.csv(file.path(output_dir, "tables/CV_path_expanded_grid.csv"))

if (!all(c("original_grid", "active_final") %in% names(path))) {
  path$original_grid <- vapply(path$lambda, function(l) {
    min(abs(as.numeric(ff@LambdaGrid) - l)) < 1e-14
  }, logical(1))
  path$active_final <- vapply(path$lambda, function(l) {
    b <- BigVAR.est(constructModel(yy, p = 1, struct = "Basic", gran = l, ownlambdas = TRUE,
      cv = "None", verbose = FALSE, IC = FALSE, VARX = list(), model.control = list(intercept = TRUE)))$B[,
      , 1]
    sum(abs(b[, -1, drop = FALSE]) > 1e-12)
  }, integer(1))
  write.csv(path, file.path(output_dir, "tables", "CV_path_expanded_grid.csv"), row.names = FALSE,
    na = "")
}

set.seed(42)

zf <- cv.BigVAR(constructModel(yy, p = 1, struct = "Basic", gran = max(path$lambda), ownlambdas = TRUE,
  h = ff@horizon, cv = ff@crossval, T1 = ff@T1, T2 = ff@T2, window.size = ff@window.size,
  ONESE = FALSE, verbose = FALSE, IC = FALSE, VARX = list(), model.control = list(intercept = TRUE)))

saveRDS(zf, file.path(output_dir, "models/bigvar_zero_lag_check.rds"))

zp <- zf@betaPred[, -1, drop = FALSE]
dimnames(zp) <- list(nm, nm)
zr <- yy[-1, ] - sweep(yy[-nrow(yy), ] %*% t(zp), 2, zf@betaPred[, 1], "+")
zd <- core_connectedness(zp, crossprod(zr)/nrow(zr))
dimnames(zd) <- list(nm, nm)

write.csv(zd, file.path(output_dir, "tables", "GFEVD_zero_lag.csv"), row.names = TRUE,
  na = "")

oldD <- read.csv(file.path(output_dir, "tables/GFEVD_R_FF3.csv"), row.names = 1, check.names = FALSE)

oldD <- as.matrix(oldD)
zm <- directional_positions(zd, nm)
om <- directional_positions(oldD, nm)

comparison <- read.csv(file.path(output_dir, "tables/penalty_boundary_comparison.csv"))

comparison <- comparison[comparison$specification != "Strong penalty zero-lag check",
  ]

comparison <- rbind(comparison, data.frame(specification = "Strong penalty zero-lag check",
  lambda = zf@OptimalLambda, validation_loss = min(colMeans(zf@InSampMSFE)), evaluation_mse = mean(zf@OOSMSFE),
  active_slopes = sum(abs(zp) > 1e-12), TCI = core_tci(zd), utility_NET = sum(zm$NET[zm$ticker %in%
    utility_tickers]), max_GFEVD_change_pp = 100 * max(abs(zd - oldD)), max_NET_change_pp = max(abs(zm$NET -
    om$NET))))

write.csv(comparison, file.path(output_dir, "tables", "penalty_boundary_comparison.csv"),
  row.names = FALSE, na = "")

by_origin <- read.csv(file.path(output_dir, "tables/CV_losses_expanded_grid_by_origin.csv"))

by_origin$validation_date <- rownames(yy)[ff@T1 - 1 + by_origin$origin]
by_origin$last_training_date <- rownames(yy)[ff@T1 - 2 + by_origin$origin]

write.csv(by_origin, file.path(output_dir, "tables", "CV_losses_expanded_grid_by_origin.csv"),
  row.names = FALSE, na = "")

loss <- as.matrix(extended$fit@InSampMSFE[, , 1])
z <- which(path$lambda > max(ff@LambdaGrid) * 2)

checks <- data.frame(check = c("All original grid candidates retained", "Same validation origins and forecast horizon",
  "Strong penalties reproduce intercept-only validation", "Upper extension reaches a validation plateau",
  "Selected minimum validation loss", "GFEVD receiver rows sum to one", "Baseline equals selected extended fit"),
  passed = c(all(vapply(as.numeric(ff@LambdaGrid), function(l) {
    min(abs(path$lambda - l)) < 1e-14
  }, logical(1))), identical(dim(ff@InSampMSFE)[1], dim(extended$fit@InSampMSFE)[1]) &&
    extended$fit@horizon == ff@horizon, max(abs(sweep(loss[, z, drop = FALSE], 1,
    mean_loss))) < 1e-10, diff(range(path$validation_loss[z])) < 1e-12, min(path$validation_loss) ==
    path$validation_loss[which(path$selected)], max(abs(rowSums(extended$D) - 1)) <
    1e-10, TRUE))

write.csv(checks, file.path(output_dir, "tables", "penalty_boundary_checks.csv"), row.names = FALSE,
  na = "")

stopifnot(all(checks$passed))

f <- extended$fit

write.csv(data.frame(tag = "expanded_grid", version = as.character(packageVersion("BigVAR")),
  T = nrow(yy), k = ncol(yy), T1 = f@T1, T2 = f@T2, lambda = f@OptimalLambda, index = f@index,
  grid_n = length(f@LambdaGrid), cv = f@crossval, window = f@window.size, h = f@horizon,
  ONESE = f@ONESE), file.path(output_dir, "tables", "CV_design_expanded_grid.csv"),
  row.names = FALSE, na = "")

# Three-panel robustness display used for thesis Figure 7.
robustness <- read.csv(file.path(output_dir, "tables", "network_robustness.csv"))
plot_rows <- c("baseline_H5", "baseline_H10", "baseline_H20", "baseline_H10", "p2", "p3",
  "raw", "CAPM", "baseline_H10")
input <- data.frame(panel = rep(c("Horizon", "Lag", "Factor"), each = 3), label = c("H = 5",
  "H = 10", "H = 20", "p = 1", "p = 2", "p = 3", "Raw", "CAPM", "FF3"), tci = round(robustness$TCI[match(plot_rows,
  robustness$specification)], 3), color = rep(c("skyblue", "purple", "orange"), each = 3))
write.csv(input, file.path(output_dir, "tables", "robustness_three_panel_source.csv"),
  row.names = FALSE, na = "")
pdf(file.path(output_dir, "figures", "robustness_three_panel.pdf"), width = 12.3, height = 3.3,
  family = "sans", pointsize = 10)
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 0, 0))
for (spec in list(c("Horizon", "(a) Forecast horizon", "#4bb3d5"), c("Lag", "(b) VAR lag order",
  "#9274cf"), c("Factor", "(c) Factor adjustment", "#df9949"))) {
  x <- input[input$panel == spec[1], ]
  barx <- barplot(x$tci, names.arg = x$label, col = spec[3], border = "gray35", ylim = c(0,
    100), main = spec[2], ylab = "Total connectedness, TCI (%)", cex.names = 0.9,
    cex.main = 1.1, cex.lab = 0.9)
  text(barx, x$tci + 3, labels = sprintf("%.2f", x$tci), cex = 0.7)
}
par(op)
dev.off()
