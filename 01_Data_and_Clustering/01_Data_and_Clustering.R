setwd("/Users/amyhu/Research/AI_Energy_Thesis_Code")
# 01 Data and Clustering
# Change the working-directory line above for another computer.
data_dir <- "09_Shared_Data/inputs"
output_dir <- "01_Data_and_Clustering/outputs"
for (folder_name in c("tables", "models", "figures", "data")) {
  dir.create(file.path(output_dir, folder_name), recursive = TRUE, showWarnings = FALSE)
}
source("10_Shared_Functions/Shared_Functions.R")

heatmap_figure_dir <- file.path(output_dir, "figures")

# Daily returns and FF3 adjustment
# Outer-join adjusted prices; forward-fill prices only; log-difference. Subtract RF,
# divide the published factors by 100, and estimate FF3 separately in each period.

stock_tickers <- sort(unlist(sector_groups, use.names = FALSE))

stopifnot(length(stock_tickers) == 94, !anyDuplicated(stock_tickers))

price_series <- lapply(stock_tickers, function(s) {
  q <- read_adjusted_prices(paste0("raw_daily_", s, ".csv"))
  names(q)[2] <- s
  q
})

names(price_series) <- stock_tickers
prices <- Reduce(function(a, b) {
  merge(a, b, by = "Date", all = TRUE)
}, price_series)

prices <- prices[order(prices$Date), ]

for (j in 2:ncol(prices)) {
  v <- prices[[j]]
  for (i in seq_along(v)) {
    if (i > 1 && is.na(v[i])) {
      v[i] <- v[i - 1]
    }
  }
  prices[[j]] <- v
}

r <- log(as.matrix(prices[-1, -1])/as.matrix(prices[-nrow(prices), -1]))
rownames(r) <- as.character(prices$Date[-1])
colnames(r) <- stock_tickers

write.csv(r, file.path(output_dir, "data", "fullsample_returns_R.csv"))

sample_bounds <- list(pre_ai_returns.csv = c("2018-01-03", "2022-11-30"), post_ai_returns.csv = c("2023-09-18",
  "2026-04-24"), ai_energy_lasso_100_daily_returns.csv = c("2023-09-15", "2026-04-24"))

for (f in names(sample_bounds)) {
  bounds <- sample_bounds[[f]]
  new <- r[rownames(r) >= bounds[1] & rownames(r) <= bounds[2], , drop = FALSE]
  new <- new[, colSums(!is.finite(new)) == 0, drop = FALSE]
  expected_firms <- if (f == "pre_ai_returns.csv") 90L else 94L
  stopifnot(ncol(new) == expected_firms, all(is.finite(new)),
    !anyDuplicated(rownames(new)), all(diff(as.Date(rownames(new))) > 0))
  write.csv(new, file.path(output_dir, "data", paste0("R_", f)))
}

unzip(file.path(data_dir, "F-F_daily.zip"), exdir = file.path(output_dir, "data/FF_source"))
lines <- readLines(file.path(output_dir, "data/FF_source/F-F_Research_Data_Factors_daily.csv"))
fl <- lines[grepl("^[0-9]{8},", lines)]

ff <- read.csv(text = paste(fl, collapse = "\n"), header = FALSE, col.names = c("date",
  "Mkt_RF", "SMB", "HML", "RF"))

ff$date <- as.Date(as.character(ff$date), "%Y%m%d")
ff[, 2:5] <- ff[, 2:5]/100

write.csv(ff, file.path(output_dir, "data", "FF3_daily_R.csv"), row.names = FALSE)

for (f in c("pre_ai_returns.csv", "post_ai_returns.csv", "ai_energy_lasso_100_daily_returns.csv")) {
  ret <- read.csv(file.path(file.path(output_dir, "data"), paste0("R_", f)), row.names = 1,
    check.names = FALSE)
  ix <- match(as.Date(rownames(ret)), ff$date)
  stopifnot(!anyNA(ix))
  x <- as.matrix(ff[ix, c("Mkt_RF", "SMB", "HML")])
  rf <- ff$RF[ix]
  e <- factor_residuals(sweep(as.matrix(ret), 1, rf), x)
  dimnames(e) <- dimnames(ret)
  write.csv(e, file.path(output_dir, "data", paste0("R_FF3_", f)))
}

sampleinfo <- do.call(rbind, lapply(c("pre_ai_returns.csv", "post_ai_returns.csv", "ai_energy_lasso_100_daily_returns.csv"),
  function(f) {
    d <- read.csv(file.path(file.path(output_dir, "data"), paste0("R_", f)), row.names = 1,
      check.names = FALSE)
    data.frame(sample = f, n = nrow(d), k = ncol(d), start = min(rownames(d)), end = max(rownames(d)))
  }))

write.csv(sampleinfo, file.path(output_dir, "tables", "sample_manifest.csv"), row.names = FALSE,
  na = "")

sp <- read_fred_series("raw_daily_sp500.csv")
sp <- sp[is.finite(sp$value), ]
sp$mkt <- c(NA, diff(log(sp$value)))
sp <- sp[, c("Date", "mkt")]
dg <- read_fred_series("raw_DGS10_daily.csv")
dg <- dg[is.finite(dg$value), ]
dg$dy10 <- c(NA, diff(dg$value))
dg <- dg[, c("Date", "dy10")]

saveRDS(list(sp = sp, dg = dg), file.path(output_dir, "models/market_controls.rds"))

# Mantegna-distance heatmaps
# Raw-return correlations; sqrt(2*(1-rho)); Ward.D2. Preserve the four distinct panels.

render_ward_heatmap <- function(csv_path, output_png, title_text) {
  df_ret <- read.csv(csv_path, row.names = 1)
  c_mat <- cor(df_ret, use = "pairwise.complete.obs")
  c_mat[c_mat > 1] <- 1
  c_mat[c_mat < -1] <- -1
  d_mantegna <- as.dist(sqrt(2 * (1 - c_mat)))
  hc <- hclust(d_mantegna, method = "ward.D2")
  dend <- as.dendrogram(hc)
  figure_id <- sub("[.]png$", "", basename(output_png))
  write.csv(c_mat, file.path(output_dir, "tables", paste0(figure_id, "_correlations.csv")))
  saveRDS(list(distance = "sqrt(2*(1-rho))", method = "ward.D2", merge = hc$merge, height = hc$height,
    order = hc$order, labels = hc$labels), file.path(output_dir, "models", paste0(figure_id,
    "_clustering.rds")))
  col_colors <- sapply(colnames(c_mat), ticker_color)
  png(output_png, width = 2400, height = 2400, res = 200)
  heatmap(c_mat, Rowv = dend, Colv = dend, scale = "none", margins = c(5, 5), col = colorRampPalette(c("blue",
    "white", "red"))(100), RowSideColors = col_colors, ColSideColors = col_colors,
    main = title_text, cexRow = 0.5, cexCol = 0.5)
  legend("topright", legend = c("AI & Tech", "Data Centers", "Nuclear & IPP", "Utilities",
    "Traditional Energy", "Renewables"), fill = c("#1f77b4", "#9467bd", "#2ca02c",
    "#ffbb78", "#d62728", "#e377c2"), bty = "n", cex = 0.8)
  dev.off()
}

render_ward_heatmap(file.path(data_dir, "pre_ai_returns.csv"), file.path(heatmap_figure_dir,
  "pre_ai_heatmap_ward.png"), "Pre-AI Regime Stock Clustering (Jan 2018 - Nov 2022, Mantegna / Ward.D2)")

render_ward_heatmap(file.path(data_dir, "post_ai_returns.csv"), file.path(heatmap_figure_dir,
  "post_ai_heatmap_ward.png"), "Post-ChatGPT Regime Stock Clustering (Sep 2023 - Apr 2026, Mantegna / Ward.D2)")

render_ward_heatmap(file.path(data_dir, "fullsample_returns.csv"), file.path(heatmap_figure_dir,
  "fullsample_90stock_heatmap_ward.png"), "Full Continuous Sample Clustering (2018 - 2026, 90 Stocks, Mantegna / Ward.D2)")

render_ward_heatmap(file.path(data_dir, "ai_energy_lasso_100_daily_returns.csv"), file.path(heatmap_figure_dir,
  "post_ai_94stock_heatmap_ward.png"), "High-Dimensional Stock Clustering (Post-ChatGPT 94 Stocks, Mantegna / Ward.D2)")

# Research definitions for descriptive thesis Tables 1 and 2.
sector_definition <- do.call(rbind, lapply(names(sector_groups), function(group) {
  data.frame(group = group, ticker = sector_groups[[group]])
}))
write.csv(sector_definition, file.path(output_dir, "tables", "Sector_Definitions.csv"),
  row.names = FALSE, na = "")
portfolio_definition <- list(chip_hardware = chip_tickers, hyperscalers = hyperscaler_tickers,
  merchant_generators = c("VST", "NRG"), regulated_utilities = utility_tickers)
portfolio_definition <- do.call(rbind, lapply(names(portfolio_definition), function(group) {
  data.frame(group = group, ticker = portfolio_definition[[group]])
}))
write.csv(portfolio_definition, file.path(output_dir, "tables", "Portfolio_Definitions.csv"),
  row.names = FALSE, na = "")
