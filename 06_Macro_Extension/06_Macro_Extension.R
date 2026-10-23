setwd("/Users/amyhu/Research/AI_Energy_Thesis_Code")
# 06 Macro Extension
# Change the working-directory line above for another computer.
data_dir <- "09_Shared_Data/inputs"
output_dir <- "06_Macro_Extension/outputs"
for (folder_name in c("tables", "models", "figures", "data")) {
  dir.create(file.path(output_dir, folder_name), recursive = TRUE, showWarnings = FALSE)
}
source("10_Shared_Functions/Shared_Functions.R")
library(glmnet)
library(sandwich)

# Rolling grid tension and monthly INDPRO predictions
# 150 days, 119 training pairs, 30 validation pairs. Window-specific FF3 and
# equation-specific chronological glmnet tuning remain distinct from full-sample BigVAR.
# Complete-month h=6 origins are May 2024–September 2025: 17 observations,
# targeting November 2024–March 2026. Do not expand the sample or change timing.

stock_returns <- as.matrix(read.csv(file.path(data_dir, "R_ai_energy_lasso_100_daily_returns.csv"),
  row.names = 1, check.names = FALSE))

archived_residuals <- as.matrix(read.csv(file.path(data_dir, "dy_ff_residuals.csv"), row.names = 1,
  check.names = FALSE))

stock_returns <- stock_returns[rownames(archived_residuals), colnames(archived_residuals)]

factor_data <- read.csv(file.path(data_dir, "FF3_daily_R.csv"))
factor_data$date <- as.Date(factor_data$date)

factor_data <- factor_data[match(as.Date(rownames(stock_returns)), factor_data$date),
  ]

factor_design <- cbind(1, as.matrix(factor_data[, c("Mkt_RF", "SMB", "HML")]))
excess <- sweep(stock_returns, 1, factor_data$RF)
k <- ncol(stock_returns)
window_length <- 150L
validation_pairs <- 30L
utility_columns <- which(colnames(stock_returns) %in% utility_tickers)

set.seed(42)

archived_penalties <- sapply(seq_len(k), function(i) {
  cv.glmnet(archived_residuals[-nrow(archived_residuals), ], archived_residuals[-1,
    i], nfolds = 5)$lambda.1se
})

write.csv(data.frame(ticker = colnames(stock_returns), lambda = archived_penalties), file.path(output_dir,
  "tables", "archived_full_sample_lambdas.csv"), row.names = FALSE, na = "")

estimate_rolling_window <- function(e, mode, window_length = 150) {
  ix <- (e - window_length + 1):e
  if (mode == "window_FF3") {
    train_days <- seq_len(window_length - validation_pairs)
    fc <- qr.coef(qr(factor_design[ix[train_days], , drop = FALSE]), excess[ix[train_days],
      , drop = FALSE])
    selection_y <- excess[ix, , drop = FALSE] - factor_design[ix, , drop = FALSE] %*%
      fc
    final_y <- factor_residuals(excess[ix, , drop = FALSE], factor_design[ix, -1,
      drop = FALSE])
  }
  else {
    selection_y <- final_y <- archived_residuals[ix, , drop = FALSE]
  }
  x <- final_y[-window_length, , drop = FALSE]
  z <- final_y[-1, , drop = FALSE]
  xs <- selection_y[-window_length, , drop = FALSE]
  zs <- selection_y[-1, , drop = FALSE]
  n <- nrow(x)
  tr <- seq_len(n - validation_pairs)
  va <- (n - validation_pairs + 1):n
  phi <- matrix(0, k, k)
  er <- matrix(0, n, k)
  lams <- numeric(k)
  for (i in seq_len(k)) {
    if (mode == "original") {
      lam <- archived_penalties[i]
    }
    else {
      path <- glmnet(xs[tr, , drop = FALSE], zs[tr, i], alpha = 1, standardize = TRUE)
      losses <- colMeans((predict(path, newx = xs[va, , drop = FALSE]) - zs[va,
        i])^2)
      lam <- path$lambda[which.min(losses)]
    }
    fit <- glmnet(x, z[, i], alpha = 1, lambda = lam, standardize = TRUE)
    phi[i, ] <- as.numeric(coef(fit))[-1]
    er[, i] <- z[, i] - as.numeric(predict(fit, newx = x))
    lams[i] <- lam
  }
  sig <- crossprod(er)/n
  dd <- core_connectedness(phi, sig)
  net <- directional_positions(dd, colnames(stock_returns))$NET
  data.frame(Date = rownames(stock_returns)[e], Util_Net_Sum = sum(net[utility_columns]),
    nonzero = sum(phi != 0), lambda_min = min(lams), lambda_median = median(lams),
    lambda_max = max(lams), first_train_date = rownames(stock_returns)[ix[1]], last_validation_date = rownames(stock_returns)[e])
}

daily_network_indices <- list()

for (mode in c("original", "within_window", "window_FF3")) {
  rows <- vector("list", nrow(stock_returns) - window_length + 1L)
  workers <- if (.Platform$OS.type == "windows") {
    1L
  }
  else {
    4L
  }
  ends <- window_length:nrow(stock_returns)
  for (start in seq(1L, length(ends), by = 25L)) {
    ids <- start:min(start + 24L, length(ends))
    batch <- parallel::mclapply(ends[ids], function(e) {
      estimate_rolling_window(e, mode)
    }, mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
    stopifnot(!any(vapply(batch, inherits, logical(1), "try-error")))
    rows[ids] <- batch
  }
  d <- do.call(rbind, rows)
  d$Date <- as.Date(d$Date)
  write.csv(d, file.path(output_dir, "tables", paste0("daily_utility_net_", mode, ".csv")),
    row.names = FALSE, na = "")
  daily_network_indices[[mode]] <- d
}

monthly <- function(d) {
  d$YearMonth <- format(d$Date, "%Y-%m")
  m <- aggregate(Util_Net_Sum ~ YearMonth, d, mean)
  m$n_days <- as.integer(table(d$YearMonth)[m$YearMonth])
  m$S_m <- as.numeric(scale(m$Util_Net_Sum))
  first <- m$YearMonth[1]
  last <- tail(m$YearMonth, 1)
  lastdate <- max(d$Date)
  nextmonth <- as.Date(paste0(add_months(last, 1), "-01"))
  lastweekday <- nextmonth - 1
  while (weekdays(lastweekday) %in% c("Saturday", "Sunday")) {
    lastweekday <- lastweekday - 1
  }
  m$complete_month <- m$YearMonth != first & (m$YearMonth != last | lastdate >= lastweekday)
  m
}

ip <- read_fred_series("raw_monthly_INDPRO.csv")
ip$YearMonth <- format(ip$Date, "%Y-%m")
ip$growth <- c(NA, diff(log(ip$value))) * 100
etf <- function(f, name) {
  d <- read_adjusted_prices(f)
  d$YearMonth <- format(d$Date, "%Y-%m")
  m <- aggregate(value ~ YearMonth, d, function(v) {
    tail(v, 1)
  })
  m[[name]] <- c(NA, diff(log(m$value))) * 100
  m[, c("YearMonth", name)]
}

xlu <- etf("raw_ai_hardware_xlu_daily.csv", "XLU")
smh <- etf("raw_ai_hardware_smh_daily.csv", "SMH")
results <- list()
idx <- 0
loo <- list()
lx <- 0
forecasts <- list()
fx <- 0

for (mode in names(daily_network_indices)) {
  m <- monthly(daily_network_indices[[mode]])
  write.csv(m, file.path(output_dir, "tables", paste0("monthly_Sm_", mode, ".csv")),
    row.names = FALSE, na = "")
  d <- Reduce(function(a, b) {
    merge(a, b, by = "YearMonth", all.x = TRUE)
  }, list(m, xlu, smh))
  d$growth <- ip$growth[match(d$YearMonth, ip$YearMonth)]
  for (completemonth in c(FALSE, TRUE)) {
    for (timing in c("same_month_ip", "lagged_ip")) {
      for (h in 1:6) {
        z <- d
        if (completemonth) {
         z <- z[z$complete_month, ]
        }
        z$y <- ip$growth[match(add_months(z$YearMonth, h), ip$YearMonth)]
        z$ip_control <- ip$growth[match(add_months(z$YearMonth, if (timing ==
         "same_month_ip") {
         0
        }
        else {
         -1
        }), ip$YearMonth)]
        z <- z[complete.cases(z[, c("y", "Util_Net_Sum", "ip_control", "XLU",
         "SMH")]), ]
        z$S_m <- as.numeric(scale(z$Util_Net_Sum))
        f <- lm(y ~ S_m + ip_control + XLU + SMH, data = z)
        se <- sqrt(diag(NeweyWest(f, lag = h, prewhite = FALSE, adjust = FALSE)))
        t <- coef(f)["S_m"]/se["S_m"]
        idx <- idx + 1
        results[[idx]] <- data.frame(index = mode, complete_months = completemonth,
         controls = timing, horizon = h, n = nobs(f), beta = coef(f)["S_m"],
         se = se["S_m"], p = 2 * pt(-abs(t), df = f$df.residual), adj_r2 = summary(f)$adj.r.squared,
         start = min(z$YearMonth), end = max(z$YearMonth))
        if (mode == "window_FF3" && completemonth && timing == "lagged_ip") {
         for (j in seq_len(nrow(z))) {
          sub <- z[-j, ]
          sub$S_m <- as.numeric(scale(sub$Util_Net_Sum))
          fm <- lm(y ~ S_m + ip_control + XLU + SMH, sub)
          sx <- sqrt(diag(NeweyWest(fm, lag = h, prewhite = FALSE)))
          lx <- lx + 1
          loo[[lx]] <- data.frame(horizon = h, omitted = z$YearMonth[j], beta = coef(fm)["S_m"],
           p = 2 * pt(-abs(coef(fm)["S_m"]/sx["S_m"]), fm$df.residual))
         }
         for (j in seq_len(nrow(z))) {
          origin <- z$YearMonth[j]
          eligible <- which(add_months(z$YearMonth, h + 1) <= origin)
          if (length(eligible) < 10) {
           next
          }
          tr <- z[eligible, ]
          va <- z[j, , drop = FALSE]
          mu <- mean(tr$Util_Net_Sum)
          ss <- sd(tr$Util_Net_Sum)
          tr$S_m <- (tr$Util_Net_Sum - mu)/ss
          va$S_m <- (va$Util_Net_Sum - mu)/ss
          aug <- lm(y ~ S_m + ip_control + XLU + SMH, tr)
          benchmark <- lm(y ~ ip_control + XLU + SMH, tr)
          ar <- lm(y ~ ip_control, tr)
          fx <- fx + 1
          forecasts[[fx]] <- data.frame(horizon = h, origin = origin, n_train = nrow(tr),
           actual = va$y, augmented = as.numeric(predict(aug, va)), financial_benchmark = as.numeric(predict(benchmark,
            va)), AR = as.numeric(predict(ar, va)), historical_mean = mean(tr$y))
         }
        }
      }
    }
  }
}

res <- do.call(rbind, results)

res$BH_q <- ave(res$p, interaction(res$index, res$complete_months, res$controls), FUN = function(v) {
  p.adjust(v, "BH")
})

write.csv(res, file.path(output_dir, "tables", "all_macro_horizons.csv"), row.names = FALSE,
  na = "")

write.csv(res[res$index == "window_FF3" & res$complete_months & res$controls == "lagged_ip",
  ], file.path(output_dir, "tables", "corrected_horizons.csv"), row.names = FALSE, na = "")

write.csv(do.call(rbind, loo), file.path(output_dir, "tables", "macro_leave_one_month_out.csv"),
  row.names = FALSE, na = "")

if (length(forecasts)) {
  fc <- do.call(rbind, forecasts)
  write.csv(fc, file.path(output_dir, "tables", "macro_expanding_forecasts.csv"), row.names = FALSE,
    na = "")
  fs <- do.call(rbind, lapply(split(fc, fc$horizon), function(d) {
    data.frame(horizon = d$horizon[1], n = nrow(d), MSE_augmented = mean((d$actual -
      d$augmented)^2), MSE_financial = mean((d$actual - d$financial_benchmark)^2),
      MSE_AR = mean((d$actual - d$AR)^2), MSE_mean = mean((d$actual - d$historical_mean)^2))
  }))
  write.csv(fs, file.path(output_dir, "tables", "macro_forecast_scores.csv"), row.names = FALSE,
    na = "")
}

cv <- do.call(rbind, lapply(window_length:nrow(archived_residuals), function(e) {
  sig <- cov(archived_residuals[(e - window_length + 1):e, ])
  dd <- core_connectedness(matrix(0, k, k), sig)
  data.frame(Date = rownames(archived_residuals)[e], Util_Net_Sum = sum(directional_positions(dd,
    colnames(archived_residuals))$NET[utility_columns]))
}))

cv$Date <- as.Date(cv$Date)

write.csv(monthly(cv), file.path(output_dir, "tables", "monthly_covariance_only_diagnostic.csv"),
  row.names = FALSE, na = "")

pdf(file.path(output_dir, "figures", "monthly_Sm.pdf"), width = 10, height = 5)
par(mar = c(5, 5, 3, 1), las = 1, cex = 0.85)

m <- monthly(daily_network_indices$window_FF3)

plot(as.Date(paste0(m$YearMonth, "-01")), m$S_m, type = "l", lwd = 2, col = "#245a81",
  xlab = "Month", ylab = "Utility NET index (standard deviations)")

points(as.Date(paste0(m$YearMonth[!m$complete_month], "-01")), m$S_m[!m$complete_month],
  pch = 1, col = "red")

abline(h = 0, lty = 3)

dev.off()

q <- res[res$index == "window_FF3" & res$complete_months & res$controls == "lagged_ip",
  ]

pdf(file.path(output_dir, "figures", "macro_horizons.pdf"), width = 8, height = 5)
par(mar = c(5, 5, 3, 1), las = 1, cex = 0.85)

plot(q$horizon, q$beta, ylim = range(q$beta - 2.4 * q$se, q$beta + 2.4 * q$se), pch = 19,
  xlab = "Forecast horizon (months)", ylab = "INDPRO monthly growth, percentage points")

cr <- qt(0.975, q$n - 5) * q$se

arrows(q$horizon, q$beta - cr, q$horizon, q$beta + cr, angle = 90, code = 3, length = 0.05)

abline(h = 0, lty = 3)

dev.off()

