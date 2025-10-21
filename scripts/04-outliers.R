# function from phiper, computing alpha diversity (richness + Shannon)
out_alpha <- compute_alpha_diversity(
  ps_merged_box_bin,
  "big_group",
  "peptide_id"
)

# ---------- SETTINGS ----------
# data.frame with alpha results per group
out_tbl <- out_alpha$big_group

# variables to use for multivariate outlier detection
vars <- c("richness")

# small regularization added to covariance when near-singular
small_eps <- 1e-8

# histogram / density bins
bins <- 30

# output directory for plots and tables
out_dir <- "results/outliers"
dir.create(out_dir, showWarnings = FALSE)

# significance for cutoff (used for mahalanobis thresholding)
alpha <- 0.01
p_vars <- length(vars)

# ---------- helper: compute mahalanobis per-group (store d2 and sqrt) ---------
# df_group: data.frame for single group
# vars: character vector of column names to use
# small_eps: tiny value to regularize covariance if needed
compute_group_maha <- function(df_group, vars, small_eps = 1e-8) {
  # copy input to avoid side-effects (we are in scope, but to be 100% safe)
  df <- df_group

  # find rows with complete data for the chosen vars
  cc_idx <- complete.cases(df[, vars])
  X_cc <- as.matrix(df[cc_idx, vars])
  p <- ncol(X_cc)

  # initialize output vectors as NA
  maha2 <- rep(NA_real_, nrow(df))
  maha <- rep(NA_real_, nrow(df))
  pval <- rep(NA_real_, nrow(df))

  # if no complete rows, return df with NA maha columns
  if (nrow(X_cc) == 0) {
    # originally a warning was emitted here; removed per request
    df$maha2 <- maha2
    df$maha <- maha
    df$maha_pval <- pval
    return(df)
  }

  # group mean and covariance (with safe error handling)
  mu <- colMeans(X_cc)
  S <- tryCatch(cov(X_cc), error = function(e) NA)

  # if covariance is invalid or singular, add small diagonal regularization
  if (any(is.na(S)) || det(S) <= .Machine$double.eps) {
    S <- S + diag(small_eps, p)
  }

  # squared mahalanobis distances for complete cases
  d2 <- mahalanobis(X_cc, center = mu, cov = S)

  # fill results back into full-length vectors
  maha2[cc_idx] <- d2
  maha[cc_idx] <- sqrt(d2)
  pval[cc_idx] <- pchisq(d2, df = p, lower.tail = FALSE)

  # attach columns and return
  df$maha2 <- maha2
  df$maha <- maha
  df$maha_pval <- pval
  return(df)
}

# --- safe file name helper ---
# replace whitespace and non-alphanumeric characters with underscores,
# limit to first 80 characters
clean_name <- function(x) {
  x2 <- as.character(x)
  x2 <- gsub("\\s+", "_", x2)
  x2 <- gsub("[^A-Za-z0-9_\\-]", "_", x2)
  substr(x2, 1, 80)
}

# ------ function: make and save plots for one group (separate files) ----------
# df_group: data.frame for a single big_group
# group_name: label used in titles and filenames
# suffix: optional suffix for filenames (e.g., "_lt4000")
make_and_save_for_group <- function(df_group, group_name, suffix = "") {
  # filter out na richness for plotting hist/qq
  df_rich <- df_group %>% filter(!is.na(richness))
  safe_g <- clean_name(group_name)

  # histogram + density for richness
  p_hist <- ggplot(df_rich, aes(x = richness)) +
    geom_histogram(aes(y = ..density..),
      bins = bins, alpha = 0.5,
      color = "black", fill = "grey70"
    ) +
    geom_density(size = 1) +
    labs(
      title = paste0(
        "Histogram + density - richness (", group_name,
        suffix, ")"
      ),
      x = "Richness", y = "Density"
    ) +
    theme_phip() +
    theme(text = element_text(size = 25))
  file_hist <- file.path(out_dir, paste0(
    "hist_richness_", safe_g,
    suffix, ".png"
  ))
  ggsave(file_hist, p_hist, width = 5, height = 5, dpi = 300, bg = "white")

  # qq-plot for richness (vs. the normal distribution, as the counts are already
  # aggregated to the subject level)
  p_qq <- ggplot(df_rich, aes(sample = richness)) +
    stat_qq(size = 0.8) +
    stat_qq_line(color = "blue") +
    labs(
      title = paste0("QQ Plot - richness (", group_name, suffix, ")"),
      x = "Theoretical quantiles", y = "Sample quantiles"
    ) +
    theme_phip() +
    theme(text = element_text(size = 25))
  file_qq <- file.path(out_dir, paste0("qq_richness_", safe_g, suffix, ".png"))
  ggsave(file_qq, p_qq, width = 5, height = 5, dpi = 300, bg = "white")

  # compute mahalanobis distances and prepare mahalanobis plot
  df_maha <- compute_group_maha(df_group, vars = vars, small_eps = small_eps)
  df_maha_plot <- df_maha %>% filter(!is.na(maha))

  # compute cutoff in squared d^2 and sqrt(d^2)
  cut_d2 <- qchisq(1 - alpha, df = p_vars)
  cut_sqrt <- sqrt(cut_d2)
  n_flagged <- sum(df_maha$maha2 > cut_d2, na.rm = TRUE)

  # only plot if there are complete rows for maha
  if (nrow(df_maha_plot) > 0) {
    p_maha <- ggplot(df_maha_plot, aes(x = maha)) +
      geom_histogram(aes(y = ..density..),
        bins = bins, alpha = 0.5,
        color = "black", fill = "grey80"
      ) +
      geom_density(size = 1) +
      geom_vline(xintercept = cut_sqrt, color = "red", size = 1) +
      annotate("text",
        x = cut_sqrt, y = Inf,
        label = paste0("alpha=", alpha, "\ncut=", round(cut_sqrt, 3)),
        vjust = 2.7, hjust = 1.15, color = "red", size = 7.5, lineheight = 0.32
      ) +
      labs(
        title = paste0("Mahalanobis distances (sqrt) - ", group_name, suffix),
        x = "Mahalanobis distance (sqrt)", y = "density"
      ) +
      theme_phip() +
      theme(text = element_text(size = 25))
    file_maha <- file.path(out_dir, paste0("maha_", safe_g, suffix, ".png"))
    ggsave(file_maha, p_maha, width = 5, height = 5, dpi = 300, bg = "white")
  } else {
    # if no complete rows exist for mahalanobis, do nothing
  }

  # numeric summary for this group / variant
  summary <- tibble::tibble(
    big_group = group_name,
    n_total = nrow(df_group),
    n_rich_nonNA = sum(!is.na(df_group$richness)),
    n_complete_for_maha = sum(complete.cases(df_group[, vars])),
    n_flagged_alpha = n_flagged,
    cut_d2 = cut_d2,
    cut_sqrt = cut_sqrt
  )

  # return list with summary and df_maha (caller will bind results)
  return(list(summary = summary, df_maha = df_maha))
}

# ---------- prepare groups safely ----------
groups <- unique(as.character(out_tbl$big_group))

# ---------- loop and save ----------
all_summaries <- list()
all_maha <- list()
for (g in groups) {
  df_g <- out_tbl %>% filter(big_group == g)
  if (nrow(df_g) == 0) {
    # skip empty groups silently
    next
  }
  safe_g <- clean_name(g)

  # generate plots and maha table for the full group
  res_full <- make_and_save_for_group(df_g, g, suffix = "")
  all_summaries[[paste0(safe_g, "_full")]] <- res_full$summary
  all_maha[[paste0(safe_g, "_full")]] <- res_full$df_maha

  # special-case: make an additional variant for mom_milk with richness < 4000
  if (g == "mom_milk") {
    df_g_f <- df_g %>% filter(!is.na(richness) & richness < 4000)
    if (nrow(df_g_f) > 0) {
      res_lt <- make_and_save_for_group(df_g_f, g, suffix = "_lt4000")
      all_summaries[[paste0(safe_g, "_lt4000")]] <- res_lt$summary
      all_maha[[paste0(safe_g, "_lt4000")]] <- res_lt$df_maha
    } else {
      # variant had no rows after filtering; skip silently
    }
  }
}

# --- thresholds used later for manual flagging ---
th_lt <- 2.576 # threshold for mom_milk_lt4000 (sqrt mahalanobis)
th_full <- 10 # threshold for mom_milk_full (sqrt mahalanobis)

# ---------- 1) mom_milk_lt4000 (maha > th_lt) ----------
hits_lt <- all_maha[["mom_milk_lt4000"]] %>%
  filter(maha > th_lt) %>%
  select(sample_id, maha, maha2, everything())

# map hits to subject/timepoint using original ps object
map_lt <- ps_merged_box_bin %>%
  filter(sample_id %in% hits_lt$sample_id) %>%
  select(sample_id, subject_id, timepoint_factor, big_group) %>%
  distinct()

# ---------- 2) mom_milk_full (maha > th_full) ----------
hits_full <- all_maha[["mom_milk_full"]] %>%
  filter(maha > th_full) %>%
  select(sample_id, maha, maha2, everything())

map_full <- ps_merged_box_bin %>%
  filter(sample_id %in% hits_full$sample_id) %>%
  select(sample_id, subject_id, timepoint_factor, big_group) %>%
  distinct()

# ---------- bind and save summaries & maha table ----------
summary_tbl <- dplyr::bind_rows(all_summaries)
write.csv(summary_tbl, file = file.path(out_dir, "summary_by_group.csv"),
          row.names = FALSE)

maha_all <- dplyr::bind_rows(all_maha, .id = "source_group_variant")
write.csv(maha_all, file = file.path(out_dir, "maha_all_groups.csv"),
          row.names = FALSE)

