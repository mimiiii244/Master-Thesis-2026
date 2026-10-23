setwd("/Users/amyhu/Research/AI_Energy_Thesis_Code")
# 08 Prediction and PJM
# Change the working-directory line above for another computer.
data_dir <- "09_Shared_Data/inputs"
output_dir <- "08_Prediction_and_PJM/outputs"
for (folder_name in c("tables", "models", "figures", "data")) {
  dir.create(file.path(output_dir, folder_name), recursive = TRUE, showWarnings = FALSE)
}
source("10_Shared_Functions/Shared_Functions.R")
library(glmnet)
library(sandwich)

# Daily prediction and PJM diagnostics
# the chronological 60/20/20 split, training-only FF3/scaling and grids.
# PJM date windows and 15,000 placebos remain exploratory timing diagnostics.

y <- as.matrix(read.csv(file.path(data_dir, "R_FF3_ai_energy_lasso_100_daily_returns.csv"),
  row.names = 1, check.names = FALSE))
n <- nrow(y) - 1
a <- floor(0.6 * n)
b <- floor(0.8 * n)
k <- ncol(y)

stock_returns <- as.matrix(read.csv(file.path(data_dir, "R_ai_energy_lasso_100_daily_returns.csv"),
  row.names = 1, check.names = FALSE))[rownames(y), colnames(y)]

factor_data <- read.csv(file.path(data_dir, "FF3_daily_R.csv"))
factor_data <- factor_data[match(rownames(y), factor_data$date), ]
F <- cbind(1, as.matrix(factor_data[, c("Mkt_RF", "SMB", "HML")]))
ex <- sweep(stock_returns, 1, factor_data$RF)
transform_training_sample <- function(end) {
  days <- seq_len(end + 1)
  beta <- qr.coef(qr(F[days, , drop = FALSE]), ex[days, , drop = FALSE])
  yy <- ex - F %*% beta
  x <- yy[-nrow(yy), ]
  target <- yy[-1, ]
  xm <- colMeans(x[1:end, ])
  xs <- apply(x[1:end, ], 2, sd) * sqrt((end - 1)/end)
  ym <- colMeans(target[1:end, ])
  ys <- apply(target[1:end, ], 2, sd) * sqrt((end - 1)/end)
  list(X = sweep(sweep(x, 2, xm), 2, xs, "/"), Y = sweep(sweep(target, 2, ym), 2, ys,
    "/"))
}

fit_prediction_model <- function(X, Y, end, ev, mix, pen, ridge = FALSE) {
  tr <- seq_len(end)
  XX <- X[tr, , drop = FALSE]
  YY <- Y[tr, , drop = FALSE]
  if (ridge) {
    xc <- scale(XX, scale = FALSE)
    yc <- scale(YY, scale = FALSE)
    co <- t(solve(crossprod(xc) + diag(pen, k), crossprod(xc, yc)))
    inter <- colMeans(YY) - as.numeric(co %*% colMeans(XX))
  }
  else {
    co <- matrix(0, k, k)
    inter <- numeric(k)
    for (i in seq_len(k)) {
      f <- glmnet(XX, YY[, i], alpha = mix, lambda = pen, standardize = FALSE, intercept = TRUE,
        thresh = 1e-08)
      be <- as.numeric(coef(f))
      inter[i] <- be[1]
      co[i, ] <- be[-1]
    }
  }
  list(co = co, inter = inter, pred = sweep(X[ev, , drop = FALSE] %*% t(co), 2, inter,
    "+"), res = YY - sweep(XX %*% t(co), 2, inter, "+"))
}

t <- transform_training_sample(a)
v <- (a + 1):b
base <- mean(sweep(t$Y[v, , drop = FALSE], 2, colMeans(t$Y[1:a, ]))^2)
grid <- list()
ix <- 0

for (mix in c(1, 0.75, 0.5, 0.25)) {
  for (pen in c(1, 0.4, 0.2, 0.1, 0.03, 0.01, 0.003)) {
    f <- fit_prediction_model(t$X, t$Y, a, v, mix, pen)
    ix <- ix + 1
    grid[[ix]] <- data.frame(model = if (mix == 1) {
      "lasso"
    }
    else {
      "elastic_net"
    }, mix = mix, penalty = pen, validation_mse = mean((t$Y[v, ] - f$pred)^2), nonzero = sum(abs(f$co) >
      1e-10))
  }
}

for (pen in c(10, 100, 1000, 10000, 1e+05, 1e+06)) {
  f <- fit_prediction_model(t$X, t$Y, a, v, 0, pen, TRUE)
  ix <- ix + 1
  grid[[ix]] <- data.frame(model = "ridge", mix = 0, penalty = pen, validation_mse = mean((t$Y[v,
    ] - f$pred)^2), nonzero = sum(abs(f$co) > 1e-10))
}

gr <- do.call(rbind, grid)
gr$validation_relative <- gr$validation_mse/base

write.csv(gr, file.path(output_dir, "tables", "regularization_validation_grid.csv"), row.names = FALSE,
  na = "")

tt <- transform_training_sample(b)
ev <- (b + 1):n
meanpred <- matrix(colMeans(tt$Y[1:b, ]), length(ev), k, byrow = TRUE)
mse0 <- mean((tt$Y[ev, ] - meanpred)^2)

out <- list(data.frame(model = "mean_only", mix = NA, penalty = NA, nonzero = 0, test_mse = mse0,
  test_relative = 1, spectral_radius = 0, TCI_H1 = NA, TCI_H10 = NA, max_FEVD_abs_change = NA))

ix <- 1

for (mix in c(1, 0.75, 0.5, 0.25, 0)) {
  g <- gr[gr$mix == mix, ]
  g <- g[which.min(g$validation_mse), ]
  f <- fit_prediction_model(tt$X, tt$Y, b, ev, mix, g$penalty, mix == 0)
  s <- cov(f$res)
  d1 <- core_connectedness(f$co, s, H = 1)
  d10 <- core_connectedness(f$co, s)
  ix <- ix + 1
  out[[ix]] <- data.frame(model = g$model, mix = mix, penalty = g$penalty, nonzero = sum(abs(f$co) >
    1e-10), test_mse = mean((tt$Y[ev, ] - f$pred)^2), test_relative = mean((tt$Y[ev,
    ] - f$pred)^2)/mse0, spectral_radius = max(Mod(eigen(f$co, only.values = TRUE)$values)),
    TCI_H1 = core_tci(d1), TCI_H10 = core_tci(d10), max_FEVD_abs_change = max(abs(d10 -
      d1)))
}

f <- fit_prediction_model(tt$X, tt$Y, b, ev, 0.5, 0.003)
s <- cov(f$res)
d1 <- core_connectedness(f$co, s, H = 1)
d10 <- core_connectedness(f$co, s)
ix <- ix + 1

out[[ix]] <- data.frame(model = "dense_EN_diagnostic", mix = 0.5, penalty = 0.003, nonzero = sum(abs(f$co) >
  1e-10), test_mse = mean((tt$Y[ev, ] - f$pred)^2), test_relative = mean((tt$Y[ev, ] -
  f$pred)^2)/mse0, spectral_radius = max(Mod(eigen(f$co, only.values = TRUE)$values)),
  TCI_H1 = core_tci(d1), TCI_H10 = core_tci(d10), max_FEVD_abs_change = max(abs(d10 -
    d1)))

tr <- do.call(rbind, out)

write.csv(tr, file.path(output_dir, "tables", "regularization_test_results.csv"), row.names = FALSE,
  na = "")

write.csv(data.frame(part = c("training", "validation", "test"), start = rownames(y)[c(2,
  a + 2, b + 2)], end = rownames(y)[c(a + 1, b + 1, n + 1)], n = c(a, b - a, n - b)),
  file.path(output_dir, "tables", "regularization_split_dates.csv"), row.names = FALSE,
  na = "")

pdf(file.path(output_dir, "figures", "regularization_test.pdf"), width = 9, height = 5)
par(mar = c(5, 5, 3, 1), las = 1, cex = 0.85)

par(mar = c(7, 5, 3, 1), las = 2)

barplot(tr$test_relative, names.arg = paste(tr$model, tr$mix), col = "#457393", ylab = "Held-out MSE relative to mean-only")

abline(h = 1, lty = 3)

dev.off()

r <- read.csv(file.path(data_dir, "fullsample_returns_R.csv"), row.names = 1, check.names = FALSE)
mc <- readRDS("09_Shared_Data/frozen_models/market_controls.rds")

d <- data.frame(Date = as.Date(rownames(r)), merchant = rowMeans(r[, c("VST", "NRG")]),
  utility = rowMeans(r[, utility_tickers]), chips = rowMeans(r[, c("NVDA", "AMD", "AVGO",
    "AMAT", "MU", "INTC", "QCOM", "TSM", "LRCX", "KLAC", "TXN", "ADI", "MRVL")]),
  oil = rowMeans(r[, c("XOM", "CVX")]))

d$difference <- d$merchant - d$utility
d <- merge(d, mc$sp, by = "Date")
d <- merge(d, mc$dg, by = "Date")
d <- d[d$Date >= as.Date("2023-09-18") & d$Date <= as.Date("2026-04-24"), ]
d <- d[complete.cases(d), ]
d <- d[order(d$Date), ]
f <- lm(difference ~ mkt + oil + dy10 + chips, d)
d$residual <- residuals(f)
dates <- as.Date(c("2024-07-30", "2025-07-22", "2025-12-17"))
events <- sapply(dates, function(dt) {
  which(d$Date >= dt)[1]
})

eventrows <- list()
ix <- 0

for (j in seq_along(dates)) {
  for (h in c(1, 2, 3, 4, 6)) {
    q <- d[events[j] + 0:(h - 1), ]
    ix <- ix + 1
    eventrows[[ix]] <- data.frame(date = dates[j], window = paste0("[0,+", h - 1,
      "]"), num_days = h, merchant_minus_utilities_cum_logret_pp = 100 * sum(q$difference),
      conditional_abnormal_cum_logret_pp = 100 * sum(q$residual))
  }
}

write.csv(do.call(rbind, eventrows), file.path(output_dir, "tables", "event_windows.csv"),
  row.names = FALSE, na = "")

eligible <- which(seq_len(nrow(d)) > 11 & seq_len(nrow(d)) < nrow(d) - 11 & sapply(seq_len(nrow(d)),
  function(j) {
    all(abs(j - events) > 12)
  }))

set.seed(18283)

perm <- list()

for (h in c(1, 3, 5)) {
  observed <- sapply(events, function(j) {
    sum(d$residual[j + 0:(h - 1)])
  })
  pool <- replicate(15000, mean(sapply(sample(eligible, 3), function(j) {
    sum(d$residual[j + 0:(h - 1)])
  })))
  ci <- quantile(pool, c(0.025, 0.975))
  perm[[as.character(h)]] <- data.frame(days = h, mean_cum_abnormal_pp = 100 * mean(observed),
    null_2_5 = 100 * ci[1], null_97_5 = 100 * ci[2], two_sided_permutation_p = mean(abs(pool) >=
      abs(mean(observed))), events = 3)
}

write.csv(do.call(rbind, perm), file.path(output_dir, "tables", "event_permutation.csv"),
  row.names = FALSE, na = "")

br <- list()

for (date in c("2024-03-01", "2024-05-01", "2024-07-30", "2024-10-01", "2025-01-02", "2025-04-01",
  "2025-07-22", "2025-10-01", "2025-12-17")) {
  q <- d
  q$post <- as.integer(q$Date >= as.Date(date))
  f <- lm(difference ~ post * (chips + mkt + oil + dy10), q)
  z <- hac_normal_inference(f)
  z <- z[z$term == "post:chips", ]
  br[[date]] <- data.frame(break_date = date, pre_n = sum(q$post == 0), post_n = sum(q$post ==
    1), chip_slope_pre = coef(f)["chips"], change = z$beta, se = z$se, p = z$p, CI_low = z$low,
    CI_high = z$high)
}

write.csv(do.call(rbind, br), file.path(output_dir, "tables", "breakpoint_sensitivity.csv"),
  row.names = FALSE, na = "")

