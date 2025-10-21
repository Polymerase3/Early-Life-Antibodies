################################################################################
########################### COMPUTING STABILITIES ##############################
################################################################################
# actually same as with the heatmaps but only for the kids; main workhorse of
# the whole script
stab_df_kids <- compute_repertoire_similarity(
  filter(ps_merged_box_bin, big_group == "kid_serum"),
  group_col = "big_group",
  time_col = "timepoint_recoded",
  similarity = "kulczynski",
  dyade_col = "dyade",
  time_mode = "pairwise",
  time_pairing = c("same", "cross"),
  subject_mode = "self",
  drop_self_time = FALSE,
  overwrite = TRUE,
  verbose = TRUE
)

# recoding the time_left and time_right back to the usual, descriptive labels
# small mapping tibble
con <- dbplyr::remote_con(stab_df_kids)
tbl_name <- dbplyr::remote_name(stab_df_kids)

# map T -> P/M
map_t2p <- tibble::tibble(
  from = c("T2", "T6", "T8"),
  to   = c("M0", "M3", "M12")
)

# temporary table with map
copy_to(con, map_t2p, name = "map_t2p_tmp", temporary = TRUE, overwrite = TRUE)

# update cols in the original stability df
DBI::dbExecute(con, sprintf(
  'UPDATE "%1$s" AS s
      SET time_left = COALESCE(m.to, s.time_left)
    FROM map_t2p_tmp AS m
   WHERE s.time_left = m.from;', tbl_name
))

DBI::dbExecute(con, sprintf(
  'UPDATE "%1$s" AS s
      SET time_right = COALESCE(m.to, s.time_right)
    FROM map_t2p_tmp AS m
   WHERE s.time_right = m.from;', tbl_name
))

DBI::dbExecute(con, sprintf('ANALYZE "%s";', tbl_name))

# ----------------------- add the metadata -------------------------------------
# add the meta:
meta_kids <- readRDS("data/meta_kids.rds")

## 0) get the DuckDB connection behind stab_df_kids
con <- dbplyr::remote_con(stab_df_kids)

## 1) stage meta_kids into DuckDB (temp)
meta_kids_db <- copy_to(
  dest = con,
  df = meta_kids,
  name = unique(DBI::SQL(paste0(
    "meta_kids_tmp_",
    as.integer(runif(1, 1e6, 9e6))
  ))),
  overwrite = TRUE,
  temporary = TRUE
)

## 2) append meta columns by subject_left (since left==right)
stab_df_kids_with_meta <- stab_df_kids %>%
  left_join(meta_kids_db, by = c("subject_left" = "subject_id"))

# done: `stab_df_kids_with_meta` is a lazy DuckDB table with meta columns
# appended

# --------------- analyze and prepare the data for plotting --------------------
## set globals for bootstrap and contrasts
set.seed(1234)
B <- 5000L
contrast_levels <- c("M0 vs M3", "M0 vs M12", "M3 vs M12")

# --- make a time rank map in DuckDB (DB-safe) ---
con <- dbplyr::remote_con(stab_df_kids_with_meta)
time_map <- tibble::tibble(time = c("M0", "M3", "M12"), rnk = c(1L, 2L, 3L))
time_map_db <- copy_to(con, time_map,
  name = paste0(
    "time_map_tmp_",
    as.integer(runif(1, 1e6, 9e6))
  ),
  temporary = TRUE, overwrite = TRUE
)

# --- build subject-level stability per ordered contrast (all DB-safe) ---
df_pairs <- stab_df_kids_with_meta %>%
  left_join(time_map_db, by = c("time_left" = "time")) %>%
  rename(r_left = rnk) %>%
  left_join(time_map_db, by = c("time_right" = "time")) %>%
  rename(r_right = rnk) %>%
  mutate(
    tmin = if_else(r_left <= r_right, time_left, time_right),
    tmax = if_else(r_left <= r_right, time_right, time_left)
  )

stab_subj_contrast <- df_pairs %>%
  group_by(subject_left, tmin, tmax) %>%
  summarise(stability = mean(similarity, na.rm = TRUE), .groups = "drop") %>%
  collect() %>%
  mutate(contrast = paste(tmin, "vs", tmax)) %>% # safe now (in R)
  select(subject_left, contrast, stability) %>%
  filter(contrast %in% contrast_levels)

# --- meta (in R) ---
meta_cols <- setdiff(colnames(meta_kids), "subject_id")
meta_subj <- meta_kids %>% rename(subject_left = subject_id)

df <- stab_subj_contrast %>%
  left_join(meta_subj, by = "subject_left") %>%
  mutate(across(where(is.character), ~ dplyr::na_if(.x, "NA")))

# ---------- helpers ----------
is_discrete <- function(x) is.character(x) || is.factor(x)
non_na_nlevels <- function(x) length(unique(x[!is.na(x)]))

meta_info <- tibble::tibble(
  var = meta_cols,
  is_num = purrr::map_lgl(var, ~ is.numeric(df[[.x]])),
  is_cat = purrr::map_lgl(var, ~ is_discrete(df[[.x]])),
  nlev = purrr::map_int(var, ~ if (is_discrete(df[[.x]])) {
    non_na_nlevels(df[[.x]])
  } else {
    NA_integer_
  })
) %>%
  mutate(kind = dplyr::case_when(
    is_cat & nlev == 2 ~ "binary",
    is_cat & nlev >= 3 ~ "multicat",
    is_num ~ "numeric",
    TRUE ~ "skip"
  ))

# different methods/functions for binary/cont/multicat
vars_binary <- meta_info %>%
  filter(kind == "binary") %>%
  pull(var)

vars_multicat <- meta_info %>%
  filter(kind == "multicat") %>%
  pull(var)
vars_multicat <- vars_multicat[-1]

vars_numeric <- meta_info %>%
  filter(kind == "numeric") %>%
  pull(var)

boot_diff_means <- function(y, g, B = 2000L, eps = 0) {
  # eps lets u add a tiny jitter if u ever want to avoid zero division;
  # default = 0 (strict)
  ok <- !(is.na(y) | is.na(g))
  y <- y[ok]
  g <- droplevels(as.factor(g[ok]))
  if (nlevels(g) != 2L) {
    return(tibble::tibble())
  }

  lv <- levels(g)
  y1 <- y[g == lv[1]]
  y2 <- y[g == lv[2]]
  n1 <- length(y1)
  n2 <- length(y2)
  if (n1 < 2 || n2 < 2) {
    return(tibble::tibble())
  }

  m1 <- mean(y1)
  m2 <- mean(y2)
  s1 <- stats::sd(y1)
  s2 <- stats::sd(y2)

  # --- formula: unweighted pooled SD
  sp <- sqrt(((s1^2) + (s2^2)) / 2)

  diff_obs <- m2 - m1
  d_obs <- if (sp > 0) diff_obs / sp else NA_real_

  # bootstrap
  bdiff <- numeric(B)
  bd <- numeric(B)
  for (b in seq_len(B)) {
    s1b <- sample(y1, replace = TRUE)
    s2b <- sample(y2, replace = TRUE)
    m1b <- mean(s1b)
    m2b <- mean(s2b)
    s1b_sd <- stats::sd(s1b)
    s2b_sd <- stats::sd(s2b)
    spb <- sqrt(((s1b_sd^2) + (s2b_sd^2)) / 2) # same formula in bootstrap

    bdiff[b] <- m2b - m1b
    bd[b] <- if (spb > 0) (m2b - m1b) / spb else NA_real_
  }

  ci_diff <- stats::quantile(bdiff, c(0.025, 0.975), names = FALSE, type = 6)
  bd_ok <- bd[is.finite(bd)]
  ci_d <- if (length(bd_ok) >= 10) {
    stats::quantile(bd_ok, c(0.025, 0.975), names = FALSE, type = 6)
  } else {
    c(NA_real_, NA_real_)
  }

  # two-sided bootstrap p for the mean difference
  p <- 2 * min(mean(bdiff >= 0), mean(bdiff <= 0))
  p <- min(max(p, 0), 1)

  tibble::tibble(
    level_a = lv[1], level_b = lv[2],
    n_a = n1, n_b = n2,
    mean_a = m1, mean_b = m2,
    sd_a = s1, sd_b = s2,
    diff = diff_obs, ci_low = ci_diff[1], ci_high = ci_diff[2],
    d = d_obs, d_ci_low = ci_d[1], d_ci_high = ci_d[2],
    p = p
  )
}

boot_pairwise_multicat <- function(y, g, B = 2000L) {
  ok <- !(is.na(y) | is.na(g))
  y <- y[ok]
  g <- droplevels(as.factor(g[ok]))
  lv <- levels(g)
  if (length(lv) < 3) {
    return(tibble::tibble())
  }
  pairs <- combn(lv, 2, simplify = FALSE)
  res <- lapply(pairs, function(pair) {
    idx <- g %in% pair
    boot_diff_means(y[idx], droplevels(g[idx]), B = B) %>%
      dplyr::mutate(level_a = pair[1], level_b = pair[2], .before = 1)
  })
  res <- Filter(function(z) nrow(z) > 0, res)
  if (!length(res)) tibble::tibble() else dplyr::bind_rows(res)
}

boot_corr <- function(y, x, B = 2000L, method = "pearson") {
  ok <- !(is.na(y) | is.na(x))
  y <- y[ok]
  x <- x[ok]
  if (length(y) < 3) {
    return(tibble::tibble())
  }
  r_obs <- suppressWarnings(stats::cor(y, x, method = method))
  if (is.na(r_obs)) {
    return(tibble::tibble())
  }
  rstar <- replicate(B, {
    idx <- sample.int(length(y), replace = TRUE)
    suppressWarnings(stats::cor(y[idx], x[idx], method = method))
  })
  ci <- stats::quantile(rstar, c(0.025, 0.975), names = FALSE, type = 6)
  p <- 2 * min(mean(rstar >= 0), mean(rstar <= 0))
  p <- min(max(p, 0), 1)
  tibble::tibble(
    n = length(y), r = r_obs, ci_low = ci[1], ci_high = ci[2],
    p = p
  )
}

# ---------- run per contrast ----------
binary_tbl <- purrr::map_dfr(contrast_levels, function(ct) {
  dct <- df %>% filter(contrast == ct)

  res <- purrr::map_dfr(vars_binary, function(v) {
    out <- boot_diff_means(dct$stability, dct[[v]], B = B)
    if (!nrow(out)) {
      tibble::tibble()
    } else {
      out %>% mutate(variable = v, contrast = ct, .before = 1)
    }
  })

  if (nrow(res) > 0) {
    res %>% mutate(p_adj = p.adjust(p, method = "BH"))
  } else {
    res
  }
})

# multicategory tests: add bh-adjusted p when applicable
multicat_tbl <- purrr::map_dfr(contrast_levels, function(ct) {
  dct <- df %>% filter(contrast == ct)

  purrr::map_dfr(vars_multicat, function(v) {
    out <- boot_pairwise_multicat(dct$stability, dct[[v]], B = B)
    if (!nrow(out)) {
      tibble::tibble()
    } else {
      out %>% mutate(variable = v, contrast = ct, .before = 1)
    }
  }) %>%
    {
      if (nrow(.) > 0 && "p" %in% names(.)) {
        mutate(., p_adj = p.adjust(.$p, method = "BH"))
      } else {
        .
      }
    }
})

# correlation tests: add bh-adjusted p when applicable
cor_tbl <- purrr::map_dfr(contrast_levels, function(ct) {
  dct <- df %>% filter(contrast == ct)

  purrr::map_dfr(vars_numeric, function(v) {
    out <- boot_corr(dct$stability, dct[[v]], B = B, method = "pearson")
    if (!nrow(out)) {
      tibble::tibble()
    } else {
      out %>% mutate(variable = v, contrast = ct, .before = 1)
    }
  }) %>%
    {
      if (nrow(.) > 0 && "p" %in% names(.)) {
        mutate(., p_adj = p.adjust(.$p, method = "BH"))
      } else {
        .
      }
    }
})

# ---------- tidy & pretty ----------
pretty_binary <- binary_tbl %>%
  dplyr::select(
    contrast, variable, level_a, level_b,
    n_a, n_b, mean_a, mean_b, sd_a, sd_b, # <-- show SDs
    diff, ci_low, ci_high,
    d, d_ci_low, d_ci_high,
    p, p_adj
  ) %>%
  dplyr::arrange(contrast, p_adj, p, variable) %>%
  dplyr::mutate(
    dplyr::across(c(
      mean_a, mean_b, sd_a, sd_b, diff, ci_low, ci_high, d,
      d_ci_low, d_ci_high
    ), ~ round(.x, 3)),
    dplyr::across(c(p, p_adj), ~ signif(.x, 3))
  )

pretty_multicat <- multicat_tbl %>%
  select(
    contrast, variable, level_a, level_b, n_a, n_b, mean_a, mean_b, diff,
    ci_low, ci_high, p, p_adj
  ) %>%
  arrange(contrast, variable, p_adj, p) %>%
  mutate(
    across(c(mean_a, mean_b, diff, ci_low, ci_high), ~ round(.x, 3)),
    across(c(p, p_adj), ~ signif(.x, 3))
  )

pretty_cor <- cor_tbl %>%
  select(contrast, variable, n, r, ci_low, ci_high, p, p_adj) %>%
  arrange(contrast, p_adj, p, variable) %>%
  mutate(
    across(c(r, ci_low, ci_high), ~ round(.x, 3)),
    across(c(p, p_adj), ~ signif(.x, 3))
  )

# results:
pretty_binary
pretty_multicat
pretty_cor

# create output dir (silent if exists) and write three tibbles to one .xlsx with
# three sheets
dir.create("results/stability_meta", recursive = TRUE, showWarnings = FALSE)

writexl::write_xlsx(
  list(
    binary   = as.data.frame(pretty_binary),
    multicat = as.data.frame(pretty_multicat),
    cor      = as.data.frame(pretty_cor)
  ),
  path = "results/stability_meta/stability_meta.xlsx"
)

################################################################################
########################### PLOTTING ###########################################
################################################################################
# defining globals
contrast_levels <- c("M0 vs M3", "M0 vs M12", "M3 vs M12")
out_dir <- "results/stability_meta"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
COLS_AB <- c("#0050a1", "#0f9c69")

# helpers + workhorse
.pstars <- function(p) {
  cut(p,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns"), right = FALSE
  )
}

plot_one_binary <- function(var_name, idx) {
  # --- levels A/B ---
  lv <- pretty_binary %>%
    filter(variable == var_name) %>%
    select(level_a, level_b) %>%
    slice(1) %>%
    tidyr::unite("lv_pair", level_a, level_b, sep = "|||") %>%
    pull(lv_pair)
  if (length(lv) == 1 && !is.na(lv)) {
    lv <- strsplit(lv, "\\|\\|\\|")[[1]]
  } else {
    lv <- df %>%
      pull(all_of(var_name)) %>%
      unique()
    lv <- lv[!is.na(lv)]
    if (length(lv) != 2) {
      return(invisible(NULL))
    }
  }
  pal_named <- stats::setNames(COLS_AB, lv)

  # --- data for the plotting ---
  d_plot <- df %>%
    select(subject_left, contrast, stability, grp = all_of(var_name)) %>%
    filter(contrast %in% contrast_levels, !is.na(grp)) %>%
    mutate(
      contrast = factor(contrast, levels = contrast_levels),
      grp = factor(grp, levels = lv)
    )
  if (nrow(d_plot) == 0L) {
    return(invisible(NULL))
  }

  # --- labels with counts ---
  n_lab <- d_plot %>%
    count(contrast, grp, name = "n") %>%
    tidyr::complete(contrast,
      grp = factor(lv, levels = lv),
      fill = list(n = 0)
    ) %>%
    tidyr::pivot_wider(names_from = grp, values_from = n) %>%
    mutate(lbl = sprintf(
      "%s\n%s: %d | %s: %d",
      as.character(contrast), lv[1], !!sym(lv[1]), lv[2],
      !!sym(lv[2])
    )) %>%
    select(contrast, lbl)
  x_labels <- setNames(n_lab$lbl, n_lab$contrast)

  # --- height of the braces ---
  r_all <- range(d_plot$stability, na.rm = TRUE)
  dy <- diff(r_all)
  if (!is.finite(dy) || dy == 0) dy <- 1
  y_pad <- 0.06 * dy
  tick_h <- 0.015 * dy
  text_off <- 0.01 * dy

  # --- 1) annotation_df in the desired format (start/end as GROUP NAMES) ---
  annotation_df <- pretty_binary %>%
    dplyr::filter(
      variable == var_name,
      contrast %in% levels(d_plot$contrast)
    ) %>%
    dplyr::distinct(contrast, p_adj) %>%
    dplyr::left_join(
      dplyr::summarise(dplyr::group_by(d_plot, contrast),
        y = max(stability, na.rm = TRUE) + y_pad,
        .groups = "drop"
      ),
      by = "contrast"
    ) %>%
    dplyr::mutate(
      start = lv[1], # group levels
      end   = lv[2],
      label = as.character(.pstars(p_adj)),
      color = var_name
    ) %>%
    dplyr::select(color, start, end, y, label, contrast)

  # --- 2) conversion of start/end -> xmin/xmax for plotting (box positions) ---
  ann_segs <- annotation_df %>%
    dplyr::mutate(
      x_ctr = as.numeric(factor(contrast, levels = levels(d_plot$contrast))),
      xmin  = x_ctr + 0.2, # left box (grp==start)
      xmax  = x_ctr - 0.2 # right box (grp==end)
    )

  # --- plot ---
  pd <- position_dodge(width = 0.8)
  pjd <- position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)

  p <- ggplot(d_plot, aes(x = contrast, y = stability, color = grp)) +
    geom_violin(aes(group = interaction(contrast, grp)),
      fill = NA, linewidth = 0.6, width = 0.9, position = pd,
      trim = TRUE,
      show.legend = FALSE
    ) +
    geom_boxplot(aes(group = interaction(contrast, grp)),
      width = 0.55, outlier.shape = NA, position = pd,
      linewidth = 0.5,
      show.legend = FALSE
    ) +
    geom_point(aes(group = grp),
      position = pjd, alpha = 0.8, size = 1.7,
      stroke = 0,
      show.legend = TRUE
    ) +
    scale_color_manual(
      values = pal_named, breaks = lv, drop = FALSE,
      name = var_name
    ) +
    scale_x_discrete(labels = x_labels) +
    labs(
      title = paste0("Stability by ", var_name), x = NULL,
      y = "Stability (Kulczynski)"
    ) +
    theme_phip() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(margin = margin(t = 6), lineheight = 0.4),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      legend.box.spacing = unit(0, "pt"),
      legend.box.margin = margin(t = -3, r = 0, b = 0, l = 0),
      legend.margin = margin(t = -3, b = 0),
      text = element_text(size = 25)
    ) +
    guides(
      color = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        override.aes = list(size = 5)
      )
    ) +
    ggsignif::geom_signif(
      data = ann_segs,
      aes(
        xmin = xmin, xmax = xmax, annotations = label, y_position = y,
        group = contrast
      ),
      textsize = 10, vjust = -0.2,
      manual = TRUE,
      inherit.aes = FALSE
    )

  out_file <- file.path("results/stability_meta", sprintf(
    "box%02d-%s.png", idx,
    var_name
  ))
  ggsave(out_file, p, width = 8, height = 5, dpi = 300, bg = "white")

  message("Saved: ", normalizePath(out_file))
  invisible(p)
}

# --- loop on all binary varaibles ---
i <- 1L
for (v in vars_binary) {
  try(plot_one_binary(v, i), silent = FALSE)
  i <- i + 1L
}
