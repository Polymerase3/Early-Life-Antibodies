## source helper functions from other file to not overcrowd this one
source("scripts/99-utils.R")

## cave! the code here is extremely repetitive and i know it - i simply didn't
## have time to deal with the issue + it works

## ------------------- WIDE MATRIX + META KIDS ---------------------------------
# converting the long data to wide matrix --> necessary to calculate the
# distances; the calculation is RAM intensive and may crash PCs with 32Gb RAM;
# i precalculated the wide matrix and cached it in the data dir, so the script
# can use the precomputed version
# path for cache
full_mat_path <- file.path("data", "full_mat.rds")

# try to load if present, if missing/corrupted - compute and save
if (file.exists(full_mat_path)) {
  full_mat <- tryCatch(
    readRDS(full_mat_path),
    error = function(e) {
      message("Cached file unreadable, recomputing: ", conditionMessage(e))
      long_to_wide_exist(ps_merged_box_bin$data_long,
        row_cols = c(
          "subject_id", "big_group",
          "timepoint_factor", "dyade"
        ),
        sparse = TRUE
      )
    }
  )

  # if we recomputed due to read error, resave
  if (!file.exists(full_mat_path) || inherits(full_mat, "dgCMatrix") &&
    length(readRDS(full_mat_path)) == 0) {
    dir.create(dirname(full_mat_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(full_mat, full_mat_path)
  }
  gc()
} else {
  # if not exists --> calculate and save
  full_mat <- long_to_wide_exist(ps_merged_box_bin$data_long,
    row_cols = c(
      "subject_id", "big_group",
      "timepoint_factor", "dyade"
    ),
    sparse = TRUE
  )
  dir.create(dirname(full_mat_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(full_mat, full_mat_path)
  gc()
}

# full_mat is now in memory and cached at tmp/full_mat.rds
# now read the previously saved meta for kids
meta_kids <- readRDS("data/meta_kids.rds")

## ------------------- KULCZYNSKI DISTANCE + META PARSER -----------------------
# ensure clean sparse binary
full_mat <- drop0(full_mat)
full_mat <- as(full_mat != 0, "dMatrix")

# calculate kulczynski between the subjects
D <- kulczynski_dist(full_mat)

# safety: make sure D is a proper distance matrix
stopifnot(identical(rownames(D), colnames(D)))

nms <- rownames(D)

# parse "subject_id__big_group__timepoint_factor__dyade"
parts <- stringr::str_split_fixed(nms, "__", 4)
parsed <- tibble(
  id_full          = nms,
  subject_id       = parts[, 1],
  big_group        = parts[, 2],
  timepoint_factor = parts[, 3],
  dyade            = parts[, 4]
)

## ------------------- tSNE: 1-3 -----------------------------------------------
# seed for reproducibility
set.seed(15092025)

# performing the tSNE
## 2-dimensional
ts2 <- Rtsne(
  D,
  dims             = 2,
  is_distance      = TRUE,
  perplexity       = 30,
  check_duplicates = FALSE,
  theta            = 0.5
)

## 3-dimensional
ts3 <- Rtsne(
  D,
  dims             = 3,
  is_distance      = TRUE,
  perplexity       = 30,
  check_duplicates = FALSE,
  theta            = 0.5
)

# extracting the group names
tsne_df_2d <- parsed %>%
  mutate(
    tSNE1 = ts2$Y[, 1],
    tSNE2 = ts2$Y[, 2]
  )

tsne_df_3d <- parsed %>%
  mutate(
    tSNE1 = ts3$Y[, 1],
    tSNE2 = ts3$Y[, 2],
    tSNE3 = ts3$Y[, 3]
  )

# paths
tsne_df_2d %<>%
  mutate(tp_num = readr::parse_number(timepoint_factor)) %>%
  arrange(dyade, tp_num)

################################################################################
############################# PLOT 1 - dyades ##################################
################################################################################
# generally, i know the code for the plotting is hella repetitive, so if
# somebody reads it, please dont judge haha; iam the only person which actually
# does the analyses rn, and unfortunately have no time to think about code
# quality or refine it, although i would really like to

pal_dyades <- randomcoloR::distinctColorPalette(
  length(unique(tsne_df_2d$dyade))
)

## without lines
ggplot(tsne_df_2d, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = factor(dyade)),
    size = 1.1,
    alpha = 0.65
  ) +
  scale_color_manual(values = pal_dyades) +
  labs(
    title = "t-SNE 2D - dyades",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2"
  ) +
  theme_phip() +
  theme(
    legend.position = "none",
    text = element_text(size = 25)
  )
ggsave("results/tsne/tsne-01.png",
  bg = "white", dpi = 300,
  height = 15, width = 15,
  unit = "cm"
)

# with lines
ggplot(tsne_df_2d, aes(tSNE1, tSNE2)) +
  geom_path(aes(group = dyade),
    linewidth = 0.25, alpha = 0.25, show.legend = FALSE
  ) +
  geom_point(aes(color = factor(dyade)),
    size = 1.1,
    alpha = 0.65
  ) +
  scale_color_manual(values = pal_dyades) +
  labs(
    title = "t-SNE 2D - dyades",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2"
  ) +
  theme_phip() +
  theme(
    legend.position = "none",
    text = element_text(size = 25)
  )
ggsave("results/tsne/tsne-01-lines.png",
  bg = "white", dpi = 300,
  height = 15, width = 15,
  unit = "cm"
)

# 3d
# 1) map T0..T8 back to descriptive labels
tp_map <- c(
  T0 = "P12", T1 = "P28", T2 = "B",  T3 = "W2",
  T4 = "M1",  T5 = "M2",  T6 = "M3", T7 = "M6", T8 = "M12"
)

tsne3 <- tsne_df_3d %>%
  mutate(
    tp_desc = dplyr::recode(timepoint_factor,
      !!!tp_map,
      .default = timepoint_factor
    ),
    tp_ord = readr::parse_number(timepoint_factor), # 0..8 for ordering
    hover_text = paste0(
      "Group: ", big_group, "<br>",
      "Timepoint: ", tp_desc, " (", timepoint_factor, ")<br>",
      "Dyade: ", dyade, "<br>",
      "Subject: ", subject_id, "<br>",
      "tSNE1: ", sprintf("%.3f", tSNE1), "<br>",
      "tSNE2: ", sprintf("%.3f", tSNE2), "<br>",
      "tSNE3: ", sprintf("%.3f", tSNE3)
    )
  )

# lines need rows ordered within each dyade; drop singleton dyades for lines
lines_df <- tsne3 %>%
  arrange(dyade, tp_ord) %>%
  group_by(dyade) %>%
  filter(n() >= 2) %>%
  ungroup()

# 2) 3D scatter without lines
p <- plot_ly(
  data = tsne3,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  marker = list(size = 6, opacity = 0.9),
  color = ~ factor(dyade),
  colors = pal_dyades,
  text = ~hover_text, hoverinfo = "text"
) %>%
  layout(
    title = "t-SNE 3D - dyades",
    scene = list(
      xaxis = list(title = "t-SNE 1"),
      yaxis = list(title = "t-SNE 2"),
      zaxis = list(title = "t-SNE 3"),
      aspectmode = "data" # keep true aspect ratio
    ),
    legend = list(title = list(text = "Dyade"))
  )
htmlwidgets::saveWidget(p, "results/tsne/tsne3d-01.html", selfcontained = TRUE)

# 2) 3D scatter with lines
# points
p <- plot_ly(
  data = tsne3,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  marker = list(size = 6, opacity = 0.9),
  color = ~ factor(dyade),
  colors = pal_dyades,
  text = ~hover_text, hoverinfo = "text",
  showlegend = FALSE
)

# lines per dyade (same colors; one trace split by dyade)
p <- add_trace(
  p,
  data = lines_df,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "lines",
  split = ~dyade, # <- this gets separate paths per dyade
  color = ~ factor(dyade), colors = pal_dyades,
  opacity = 0.35,
  hoverinfo = "none",
  showlegend = FALSE
)

p <- layout(
  p,
  title = "t-SNE 3D - dyades with paths",
  scene = list(
    xaxis = list(title = "t-SNE 1"),
    yaxis = list(title = "t-SNE 2"),
    zaxis = list(title = "t-SNE 3"),
    aspectmode = "data"
  ),
  showlegend = FALSE
)
htmlwidgets::saveWidget(p, "results/tsne/tsne3d-01-lines.html",
  selfcontained = TRUE
)

# thomas wanted to check how the tsnes look like with the lines --> but imho the
# clean version with only the colors is a lot better, the lines dont really add
# any value

################################################################################
############################# PLOT 2 - groups ##################################
################################################################################
## 2D without lines
ggplot(tsne_df_2d, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = factor(big_group)),
    size = 2.2,
    alpha = 0.65
  ) +
  scale_color_brewer(palette = "Set1", name = "Group") +
  labs(
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2"
  ) +
  theme_phip() +
  theme(
    legend.position = "top",
    text = element_text(size = 30)
  )
ggsave("results/tsne/tsne-02.png",
  bg = "white", dpi = 300, height = 10,
  width = 10, unit = "cm"
)

# 3D without lines
p <- plot_ly(
  data = tsne3,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  marker = list(size = 6, opacity = 0.9),
  color = ~ factor(big_group),
  text = ~hover_text, hoverinfo = "text"
) %>%
  layout(
    title = "t-SNE 3D - groups",
    scene = list(
      xaxis = list(title = "t-SNE 1"),
      yaxis = list(title = "t-SNE 2"),
      zaxis = list(title = "t-SNE 3"),
      aspectmode = "data" # keep true aspect ratio
    ),
    legend = list(title = list(text = "Group"))
  )
p
htmlwidgets::saveWidget(p, "results/tsne/tsne3d-02.html", selfcontained = TRUE)

################################################################################
############################# PLOT 3 - groups x time ###########################
################################################################################
# --- colors ---
combo_pal <- c(
  mom_serum__T0 = "#F6A1C7",
  mom_serum__T1 = "#E32877",
  mom_serum__T2 = "#8C0303",
  kid_serum__T2 = "#B7D68F",
  kid_serum__T6 = "#91BE8C",
  kid_serum__T8 = "#0F9C69",
  mom_milk__T3 = "#73ACF1",
  mom_milk__T4 = "#6CA5EB",
  mom_milk__T5 = "#5C97DE",
  mom_milk__T6 = "#4C8AD1",
  mom_milk__T7 = "#0063AB",
  mom_milk__T8 = "#0050A1"
)

# recode legend labels: group -> short abbr and time -> T-code, join with underscore
time_map <- c(
  T0  = "P12", T1  = "P28", T2  = "B",  T3  = "W2",
  T4  = "M1",  T5  = "M2",  T6  = "M3", T7  = "M6", T8  = "M12"
)

group_map <- c(kid_serum = "I", mom_milk = "BM", mom_serum = "M")

# labeller: input keys like "mom_serum__T3" or "mom_serum__M3" or "mom_serum__P12"
# returns "M_T3" / "BM_T6" etc. uses underscore as separator
lab_combo <- function(keys) {
  vapply(keys, function(k) {
    parts <- strsplit(k, "__", fixed = TRUE)[[1]]
    grp <- parts[1]
    tp_raw <- parts[2]

    # group abbreviation (fallback to original group name if unknown)
    grp_abbr <- if (!is.na(group_map[grp])) group_map[grp] else grp

    # time: if it's a readable label (P12/P28/...), map to T-code;
    # if it's already a T-code (e.g. "T3"), keep it; otherwise fallback to raw
    if (tp_raw %in% names(time_map)) {
      tp_code <- time_map[tp_raw]
    } else if (tp_raw %in% time_map) { # already a T-code present as value
      tp_code <- tp_raw
    } else {
      tp_code <- tp_raw
    }

    paste0(grp_abbr, "_", tp_code)
  }, character(1))
}

# build plot, using the new labeller for the legend labels
df2 <- tsne_df_2d |>
  mutate(combo = paste(big_group, timepoint_factor, sep = "__"))

ggplot(df2, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = combo), size = 0.8, alpha = 0.8) +
  scale_color_manual(
    breaks = names(combo_pal),
    labels = lab_combo(names(combo_pal)),
    values = combo_pal,
    name   = "Group and Time"
  ) +
  labs(x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_phip() +
  theme(
    legend.position = "top",
    text = element_text(size = 37),
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 35),
    legend.title.position = "top",
    legend.box.just = "left",
    legend.title.align = 0.5
  )

ggsave("results/tsne/tsne-03.png",
  bg = "white", dpi = 300,
  height = 15, width = 15, units = "cm"
)

# identify siblings
df2 %>%
  filter(subject_id %in% c("150788", "120677"))

df2 %>%
  filter(subject_id %in% c("202915", "202903"))

df2 %>%
  filter(subject_id %in% c("205217", "204122"))

# ids grouped into pairs (each vector = one pair)
pairs_list <- list(
  c("150788", "120677"), # twins
  c("202915", "202903"), # twins
  c("205217", "204122"), # twins
  c("304062", "180899"), # siblings
  c("201170", "100333"), # siblings
  c("005618", "005618") # siblings
)

# expand to a tibble: subject_id -> pair_id
highlight_df <- purrr::imap_dfr(pairs_list, ~ tibble::tibble(
  subject_id = .x,
  pair_id = paste0("pair_", .y)
))

# keep only rows present in df2 (defensive)
highlight_df <- df2 %>%
  dplyr::inner_join(highlight_df, by = "subject_id")

# choose one color per pair (tweak as you like)
pair_cols <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628")
names(pair_cols) <- paste0("pair_", seq_along(pair_cols))

# plotting: base points + highlighted pair points (different color per pair)
ggplot(df2, aes(tSNE1, tSNE2)) +
  # base points colored by combo
  geom_point(aes(color = combo), size = 0.8, alpha = 0.8) +
  # highlighted pair points on top: filled circles, colored by pair_id
  geom_point(
    data = highlight_df,
    aes(x = tSNE1, y = tSNE2, fill = pair_id),
    shape = 21, colour = "black", stroke = 0.4,
    size = 4.5, show.legend = TRUE
  ) +
  scale_fill_manual(name = "highlight pairs", values = pair_cols) +
  scale_color_manual(
    breaks = names(combo_pal),
    labels = lab_combo(names(combo_pal)),
    values = combo_pal,
    name   = "Group and Time"
  ) +
  labs(x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_phip() +
  theme(
    legend.position = "top",
    text = element_text(size = 37),
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 35),
    legend.title.position = "top",
    legend.box.just = "left",
    legend.title.align = 0.5
  )

# optional: add text labels for highlighted points (requires ggrepel)
# Uncomment if you want text labels
# + ggrepel::geom_text_repel(
#     data = highlight_df,
#     aes(tSNE1, tSNE2, label = subject_id),
#     color = "black", size = 3.5, segment.size = 0.3, show.legend = FALSE
#   )


# =========================
# 3D PLOT (fixed)
# =========================

# build a colors table directly from combo_pal
combo_colors_tbl <- tibble::tibble(
  combo   = names(combo_pal),
  col_hex = unname(combo_pal)
)

# preparing data nad ensure required columns exist
tsne31 <- tsne3 %>%
  mutate(combo = paste(big_group, timepoint_factor, sep = "__")) %>%
  left_join(combo_colors_tbl, by = "combo") %>%
  mutate(
    # Human-readable label like "mom_serum - T0 (-7 m)"
    combo_label = paste0(
      big_group, " - ",
      ifelse(timepoint_factor %in% names(tp_label_map),
        tp_label_map[timepoint_factor],
        timepoint_factor
      )
    )
  )

# legend levels and colors strictly aligned with combo_pal order
legend_lvls <- vapply(names(combo_pal), function(k) {
  p <- strsplit(k, "__", fixed = TRUE)[[1]]
  grp <- p[1]
  tp <- p[2]
  paste0(grp, " - ", ifelse(tp %in% names(tp_label_map), tp_label_map[tp], tp))
}, character(1))

legend_cols <- unname(combo_pal)

tsne31$combo_label <- factor(tsne31$combo_label, levels = legend_lvls)

# robust hover text (falls back if columns missing)
safe_first <- function(x) {
  if (is.null(x) || all(is.na(x))) {
    NA_character_
  } else {
    as.character(x)
  }
}

tsne31 <- tsne31 %>%
  mutate(
    subject_id = safe_first(if ("subject_id" %in% names(.)) subject_id else NA),
    hover_text = paste0(
      combo_label,
      ifelse(is.na(subject_id), "", paste0("<br>subject: ", subject_id)),
      "<br>x: ", signif(tSNE1, 4),
      "<br>y: ", signif(tSNE2, 4),
      "<br>z: ", signif(tSNE3, 4)
    )
  )

p <- plot_ly(
  tsne31,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~combo_label, # factor -> legend
  colors = legend_cols, # aligned to levels above
  marker = list(size = 4, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text",
  showlegend = TRUE
) %>%
  layout(
    scene = list(
      xaxis = list(title = "t-SNE 1"),
      yaxis = list(title = "t-SNE 2"),
      zaxis = list(title = "t-SNE 3"),
      aspectmode = "data"
    ),
    legend = list(title = list(text = "tSNE 3D - Group x Time"))
  )

htmlwidgets::saveWidget(p, "results/tsne/tsne3d-03.html", selfcontained = TRUE)


################################################################################
############################# PLOT 4 - time kids only ##########################
################################################################################
# thomas wanted to take a look, how the kids looked like on the tsnes with
# respect to different covariates --> i subsetted the big kulczynski distance
# matrix and then plotted the kids only aghainst the previously defined metadata
# later we will also do separate analyses only for month 12, but tbh it showed
# nothing interesting

# subsetting the distance matrix only to kids
# safety: make sure D is a proper distance matrix
stopifnot(identical(rownames(D), colnames(D))) # must be square

nms <- rownames(D)
parts <- do.call(rbind, strsplit(nms, "__", fixed = TRUE))
colnames(parts) <- c("subject_id", "big_group", "timepoint_factor", "dyade")
meta <- as.data.frame(parts, stringsAsFactors = FALSE)

idx <- meta$big_group == "kid_serum"
D_kid <- D[idx, idx, drop = FALSE]

# parse "subject_id__big_group__timepoint_factor__dyade"
parts_kid <- stringr::str_split_fixed(nms, "__", 4)
parsed_kid <- tibble(
  id_full          = nms[idx],
  subject_id       = parts[idx, 1],
  big_group        = parts[idx, 2],
  timepoint_factor = parts[idx, 3],
  dyade            = parts[idx, 4]
)

## ------------------- tSNE ------------------------------------
# seed for reproducibility
set.seed(15092025)

# performing the tSNE
## 2-dimensional
ts2_kid <- Rtsne(
  D_kid,
  dims             = 2,
  is_distance      = TRUE,
  perplexity       = 30,
  check_duplicates = FALSE,
  theta            = 0.5
)

## 3-dimensional
ts3_kid <- Rtsne(
  D_kid,
  dims             = 3,
  is_distance      = TRUE,
  perplexity       = 30,
  check_duplicates = FALSE,
  theta            = 0.5
)

# extracting the group names
tsne_df_2d_kid <- parsed_kid %>%
  mutate(
    tSNE1 = ts2_kid$Y[, 1],
    tSNE2 = ts2_kid$Y[, 2]
  ) %>%
  left_join(select(meta_kids, -dyade), by = "subject_id")

tsne_df_3d_kid <- parsed_kid %>%
  mutate(
    tSNE1 = ts3_kid$Y[, 1],
    tSNE2 = ts3_kid$Y[, 2],
    tSNE3 = ts3_kid$Y[, 3]
  ) %>%
  left_join(select(meta_kids, -dyade), by = "subject_id")

## 2d plot -------
# legend label map
tp_map <- c(
  T0 = "P12", T1 = "P28", T2 = "B", T3 = "W2",
  T4 = "M1", T5 = "M2", T6 = "M3", T7 = "M6", T8 = "M12"
)

tp_levels <- names(tp_map)

# color palette (9 levels) mapped to T0..T8
pal_tp <- setNames(brewer.pal(length(tp_levels), "Set1"), tp_levels)

df <- tsne_df_2d_kid %>%
  mutate(tp = factor(timepoint_factor, levels = tp_levels))

# show only timepoints that are present in the data
present <- levels(droplevels(df$tp))

ggplot(df, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = tp), size = 2, alpha = 1) +
  scale_color_manual(
    values = pal_tp[present],
    breaks = present,
    labels = unname(tp_map[present]),
    name   = "Timepoint"
  ) +
  labs(
    title = "t-SNE 2D - kids (colored by timepoint)",
    x = "t-SNE Dimension 1", y = "t-SNE Dimension 2"
  ) +
  theme_phip() +
  theme(
    legend.position = "top",
    text = element_text(size = 25)
  )
ggsave("results/tsne/tsne-04.png",
  bg = "white", dpi = 300, height = 15,
  width = 15, unit = "cm"
)

## 3d plot --------------
# map codes -> labels (legend text)
tp_map <- c(
  T0 = "P12", T1 = "P28", T2 = "B", T3 = "W2", T4 = "M1", T5 = "M2", T6 = "M3",
  T7 = "M6", T8 = "M12"
)
tp_levels <- names(tp_map)

# pick a palette for up to 9 levels (same as your 2D)
pal_tp <- setNames(brewer.pal(length(tp_levels), "Set1"), tp_levels)

# filter to kids and prepare factor/labels
tsne3_kid <- tsne31 %>%
  filter(big_group == "kid_serum") %>%
  mutate(tp = factor(timepoint_factor, levels = tp_levels))

present <- levels(droplevels(tsne3_kid$tp)) # T* present in data
legend_cols <- unname(pal_tp[present]) # colors in T* order
legend_lvls <- unname(tp_map[present]) # labels in same order

tsne3_kid <- tsne3_kid %>%
  mutate(tp_label = factor(tp_map[as.character(tp)], levels = legend_lvls))

# 3D scatter with legend
p <- plot_ly(
  tsne3_kid,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~tp_label, # legend uses mapped labels (P12, M1, ...)
  colors = legend_cols, # aligned to levels above
  marker = list(size = 6, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text",
  showlegend = TRUE
) %>%
  layout(
    scene = list(
      xaxis = list(title = "t-SNE 1"),
      yaxis = list(title = "t-SNE 2"),
      zaxis = list(title = "t-SNE 3"),
      aspectmode = "data"
    ),
    legend = list(title = list(text = "t-SNE 3D - kids (colored by timepoint)"))
  )

# save
htmlwidgets::saveWidget(p, "results/tsne/tsne3d-04.html", selfcontained = TRUE)

################################################################################
############################# PLOT 5 - time kids only ##########################
################################################################################
# recode + order + keep an explicit NA level
df_deliv <- df %>%
  mutate(
    delivery_label = recode(delivery_mode,
      VG = "Vaginal",
      CS = "C-section",
      .default = NA_character_
    ),
    delivery_label = fct_explicit_na(delivery_label, na_level = "Unknown"),
    delivery_label = fct_relevel(
      delivery_label, "Vaginal", "C-section",
      "Unknown"
    )
  )

# color-blind palette
pal_delivery <- c(
  "Vaginal"   = "#0072B2",
  "C-section" = "#D55E00",
  "Unknown"   = "grey70"
)

present <- levels(droplevels(df_deliv$delivery_label))

ggplot(df_deliv, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = delivery_label), size = 2, alpha = 0.9) +
  scale_color_manual(
    values = pal_delivery[present],
    breaks = present,
    name = "Delivery mode"
  ) +
  labs(
    title = "t-SNE 2D - kids (colored by delivery mode)",
    x = "t-SNE Dimension 1", y = "t-SNE Dimension 2"
  ) +
  theme_phip() +
  theme(
    legend.position = "top",
    text = element_text(size = 25)
  )
ggsave("results/tsne/tsne-05.png",
  bg = "white", dpi = 300, height = 15,
  width = 15, unit = "cm"
)

# 3d
# start from your 3D kids data (with tSNE1/2/3 + metadata)
df3 <- tsne_df_3d_kid

# recode delivery mode -> legend labels
df3 <- df3 %>%
  mutate(
    delivery_label = recode(delivery_mode,
      VG = "Vaginal",
      CS = "C-section",
      .default = NA_character_
    ),
    delivery_label = fct_explicit_na(delivery_label, na_level = "Unknown"),
    delivery_label = fct_relevel(delivery_label, "Vaginal", "C-section", "Unknown"),
    hover_text = paste0(
      "Subject: ", subject_id, "<br>",
      "Dyade: ", dyade, "<br>",
      "Delivery: ", delivery_label, "<br>",
      "Timepoint: ", timepoint_factor, "<br>",
      "tSNE1: ", sprintf("%.3f", tSNE1), "<br>",
      "tSNE2: ", sprintf("%.3f", tSNE2), "<br>",
      "tSNE3: ", sprintf("%.3f", tSNE3)
    )
  )

# color-blind–safe palette (order must match factor levels actually present)
pal_delivery <- c(
  "Vaginal" = "#0072B2", "C-section" = "#D55E00",
  "Unknown" = "grey70"
)
present <- levels(droplevels(df3$delivery_label))
cols <- unname(pal_delivery[present])

# 3D scatter with legend
p <- plot_ly(
  df3,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~delivery_label, colors = cols, # discrete mapping -> legend
  marker = list(size = 6, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text",
  showlegend = TRUE
) %>%
  layout(
    title = "t-SNE 3D - kids (by delivery mode)",
    scene = list(
      xaxis = list(title = "t-SNE 1"),
      yaxis = list(title = "t-SNE 2"),
      zaxis = list(title = "t-SNE 3"),
      aspectmode = "data"
    ),
    legend = list(title = list(text = "Delivery mode"))
  )

# save for slides/web
htmlwidgets::saveWidget(p, "results/tsne/tsne3d-05.html", selfcontained = TRUE)

################################################################################
############################# PLOT 6 - feedmode_birth kids only ################
################################################################################
df_feed <- df %>%
  mutate(
    feed_label = recode(feedmode_birth,
      "BF" = "Breastfeeding",
      "FF" = "Formula",
      "FF/BF" = "Mixed",
      .default = NA_character_
    ),
    feed_label = fct_explicit_na(feed_label, na_level = "Unknown"),
    feed_label = fct_relevel(
      feed_label, "Breastfeeding", "Mixed", "Formula",
      "Unknown"
    )
  )

# color-blind palette
pal_feed <- c(
  "Breastfeeding" = "#009E73",
  "Mixed"         = "#E69F00",
  "Formula"       = "#CC79A7",
  "Unknown"       = "grey70"
)

present <- levels(droplevels(df_feed$feed_label))

p2d <- ggplot(df_feed, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = feed_label), size = 2, alpha = 1) +
  scale_color_manual(
    values = pal_feed[present],
    breaks = present,
    name = "Feeding at birth"
  ) +
  labs(
    title = "t-SNE 2D - kids (feeding at birth)",
    x = "t-SNE Dimension 1", y = "t-SNE Dimension 2"
  ) +
  theme_phip() +
  theme(
    legend.position = "top",
    text = element_text(size = 25)
  )

p2d
ggsave("results/tsne/tsne-06.png",
  plot = p2d,
  bg = "white", dpi = 300, height = 15, width = 15, units = "cm"
)


# --- 3D plotly: kids colored by feedmode_birth ---
df3 <- tsne_df_3d_kid %>%
  mutate(
    feed_label = recode(feedmode_birth,
      "BF" = "Breastfeeding",
      "FF" = "Formula",
      "FF/BF" = "Mixed",
      .default = NA_character_
    ),
    feed_label = fct_explicit_na(feed_label, na_level = "Unknown"),
    feed_label = fct_relevel(
      feed_label, "Breastfeeding", "Mixed", "Formula",
      "Unknown"
    ),
    hover_text = paste0(
      "Subject: ", subject_id, "<br>",
      "Dyade: ", dyade, "<br>",
      "Feeding at birth: ", feed_label, "<br>",
      "Timepoint: ", timepoint_factor, "<br>",
      "tSNE1: ", sprintf("%.3f", tSNE1), "<br>",
      "tSNE2: ", sprintf("%.3f", tSNE2), "<br>",
      "tSNE3: ", sprintf("%.3f", tSNE3)
    )
  )

present3 <- levels(droplevels(df3$feed_label))
cols3 <- unname(pal_feed[present3])

p3d <- plot_ly(
  df3,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~feed_label, colors = cols3,
  marker = list(size = 6, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text",
  showlegend = TRUE
) %>%
  layout(
    title = "t-SNE 3D - kids (feeding at birth)",
    scene = list(
      xaxis = list(title = "t-SNE 1"),
      yaxis = list(title = "t-SNE 2"),
      zaxis = list(title = "t-SNE 3"),
      aspectmode = "data"
    ),
    legend = list(title = list(text = "Feeding at birth"))
  )

p3d
htmlwidgets::saveWidget(p3d, "results/tsne/tsne3d-06.html",
  selfcontained = TRUE
)

################################################################################
############################# PLOT 7 - sex kids only ###########################
################################################################################
df_sex <- df %>%
  mutate(
    sex_label = recode(infant_sex,
      male = "Male",
      female = "Female",
      .default = NA_character_
    ),
    sex_label = fct_explicit_na(sex_label, na_level = "Unknown"),
    sex_label = fct_relevel(sex_label, "Female", "Male", "Unknown")
  )

# color-blind palette
pal_sex <- c(
  "Female"  = "#CC79A7",
  "Male"    = "#56B4E9",
  "Unknown" = "grey70"
)

present <- levels(droplevels(df_sex$sex_label))

p2d <- ggplot(df_sex, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = sex_label), size = 2, alpha = 0.9) +
  scale_color_manual(
    values = pal_sex[present],
    breaks = present,
    name = "Infant sex"
  ) +
  labs(
    title = "t-SNE 2D - kids (by infant sex)",
    x = "t-SNE Dimension 1", y = "t-SNE Dimension 2"
  ) +
  theme_phip() +
  theme(
    legend.position = "top",
    text = element_text(size = 25)
  )

p2d
ggsave("results/tsne/tsne-07.png",
  plot = p2d, bg = "white", dpi = 300,
  height = 15, width = 15, units = "cm"
)


# --- 3D plotly: kids colored by infant_sex ---
df3 <- tsne_df_3d_kid %>%
  mutate(
    sex_label = recode(infant_sex,
      male = "Male",
      female = "Female",
      .default = NA_character_
    ),
    sex_label = fct_explicit_na(sex_label, na_level = "Unknown"),
    sex_label = fct_relevel(sex_label, "Female", "Male", "Unknown"),
    hover_text = paste0(
      "Subject: ", subject_id, "<br>",
      "Dyade: ", dyade, "<br>",
      "Infant sex: ", sex_label, "<br>",
      "Timepoint: ", timepoint_factor, "<br>",
      "tSNE1: ", sprintf("%.3f", tSNE1), "<br>",
      "tSNE2: ", sprintf("%.3f", tSNE2), "<br>",
      "tSNE3: ", sprintf("%.3f", tSNE3)
    )
  )

present3 <- levels(droplevels(df3$sex_label))
cols3 <- unname(pal_sex[present3])

p3d <- plot_ly(
  df3,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~sex_label, colors = cols3, # discrete mapping -> legend
  marker = list(size = 6, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text",
  showlegend = TRUE
) %>%
  layout(
    title = "t-SNE 3D - kids (by infant sex)",
    scene = list(
      xaxis = list(title = "t-SNE 1"),
      yaxis = list(title = "t-SNE 2"),
      zaxis = list(title = "t-SNE 3"),
      aspectmode = "data"
    ),
    legend = list(title = list(text = "Infant sex"))
  )

p3d
htmlwidgets::saveWidget(p3d, "results/tsne/tsne3d-07.html",
  selfcontained = TRUE
)

################################################################################
############################# PLOT 8 - feeding mode 3 kids only ################
################################################################################
df_m3 <- df %>%
  mutate(
    m3_label = case_when(
      feedmode_m3 == "BF" ~ "Breastfeeding",
      feedmode_m3 == "FF" ~ "Formula",
      feedmode_m3 == "FF/BF" ~ "Mixed",
      grepl("^Otherwise", feedmode_m3 %||% "") ~ "Other",
      TRUE ~ NA_character_
    ),
    m3_label = fct_explicit_na(m3_label, na_level = "Unknown"),
    m3_label = fct_relevel(
      m3_label, "Breastfeeding", "Mixed", "Formula",
      "Other", "Unknown"
    )
  )

# color-blind palette
pal_m3 <- c(
  "Breastfeeding" = "#009E73",
  "Mixed"         = "#E69F00",
  "Formula"       = "#CC79A7",
  "Other"         = "#9467BD",
  "Unknown"       = "grey70"
)

present <- levels(droplevels(df_m3$m3_label))

p2d <- ggplot(df_m3, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = m3_label), size = 2, alpha = 0.9) +
  scale_color_manual(
    values = pal_m3[present], breaks = present,
    name = "Feeding at 3 mo"
  ) +
  labs(
    title = "t-SNE 2D - kids (feeding at 3 months)",
    x = "t-SNE Dimension 1", y = "t-SNE Dimension 2"
  ) +
  theme_phip() +
  theme(legend.position = "top", text = element_text(size = 25))

p2d
ggsave("results/tsne/tsne-08.png",
  plot = p2d, bg = "white", dpi = 300, height = 15, width = 15,
  units = "cm"
)

# --- 3D plotly: kids colored by feedmode_m3 ---
df3 <- tsne_df_3d_kid %>%
  mutate(
    m3_label = case_when(
      feedmode_m3 == "BF" ~ "Breastfeeding",
      feedmode_m3 == "FF" ~ "Formula",
      feedmode_m3 == "FF/BF" ~ "Mixed",
      grepl("^Otherwise", feedmode_m3 %||% "") ~ "Other",
      TRUE ~ NA_character_
    ),
    m3_label = fct_explicit_na(m3_label, na_level = "Unknown"),
    m3_label = fct_relevel(
      m3_label, "Breastfeeding", "Mixed", "Formula",
      "Other", "Unknown"
    ),
    hover_text = paste0(
      "Subject: ", subject_id, "<br>",
      "Dyade: ", dyade, "<br>",
      "Feeding @3m: ", m3_label, "<br>",
      "Timepoint: ", timepoint_factor, "<br>",
      "tSNE1: ", sprintf("%.3f", tSNE1), "<br>",
      "tSNE2: ", sprintf("%.3f", tSNE2), "<br>",
      "tSNE3: ", sprintf("%.3f", tSNE3)
    )
  )

present3 <- levels(droplevels(df3$m3_label))
cols3 <- unname(pal_m3[present3])

p3d <- plot_ly(
  df3,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~m3_label, colors = cols3,
  marker = list(size = 6, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text",
  showlegend = TRUE
) %>%
  layout(
    title = "t-SNE 3D - kids (feeding at 3 months)",
    scene = list(
      xaxis = list(title = "t-SNE 1"),
      yaxis = list(title = "t-SNE 2"),
      zaxis = list(title = "t-SNE 3"),
      aspectmode = "data"
    ),
    legend = list(title = list(text = "Feeding at 3 mo"))
  )

p3d
htmlwidgets::saveWidget(p3d, "results/tsne/tsne3d-08.html",
  selfcontained = TRUE
)


################################################################################
############################# PLOT 9 - smoking kids only #######################
################################################################################

df_smoke <- df %>%
  mutate(
    smoke_label = case_when(
      is.na(smoking) ~ NA_character_,
      tolower(smoking) %in% c("yes", "y", "1", "true") ~ "Yes",
      tolower(smoking) %in% c("no", "n", "0", "false") ~ "No",
      TRUE ~ "Other"
    ),
    smoke_label = fct_explicit_na(smoke_label, na_level = "Unknown"),
    smoke_label = fct_relevel(smoke_label, "No", "Yes", "Unknown", "Other")
  )

# color-blind palette
pal_smoke <- c(
  "No"      = "#009E73",
  "Yes"     = "#D55E00",
  "Unknown" = "grey70",
  "Other"   = "#9467BD"
)

present <- levels(droplevels(df_smoke$smoke_label))

p2d <- ggplot(df_smoke, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = smoke_label), size = 2, alpha = 0.9) +
  scale_color_manual(
    values = pal_smoke[present], breaks = present,
    name = "Smoking"
  ) +
  labs(
    title = "t-SNE 2D - kids (smoking)",
    x = "t-SNE Dimension 1", y = "t-SNE Dimension 2"
  ) +
  theme_phip() +
  theme(legend.position = "top", text = element_text(size = 25))

p2d
ggsave("results/tsne/tsne-09.png",
  plot = p2d, bg = "white", dpi = 300, height = 15, width = 15,
  units = "cm"
)

# --- 3D plotly: kids colored by smoking ---
df3 <- tsne_df_3d_kid %>%
  mutate(
    smoke_label = case_when(
      is.na(smoking) ~ NA_character_,
      tolower(smoking) %in% c("yes", "y", "1", "true") ~ "Yes",
      tolower(smoking) %in% c("no", "n", "0", "false") ~ "No",
      TRUE ~ "Other"
    ),
    smoke_label = fct_explicit_na(smoke_label, na_level = "Unknown"),
    smoke_label = fct_relevel(smoke_label, "No", "Yes", "Unknown", "Other"),
    hover_text = paste0(
      "Subject: ", subject_id, "<br>",
      "Dyade: ", dyade, "<br>",
      "Smoking: ", smoke_label, "<br>",
      "Timepoint: ", timepoint_factor, "<br>",
      "tSNE1: ", sprintf("%.3f", tSNE1), "<br>",
      "tSNE2: ", sprintf("%.3f", tSNE2), "<br>",
      "tSNE3: ", sprintf("%.3f", tSNE3)
    )
  )

present3 <- levels(droplevels(df3$smoke_label))
cols3 <- unname(pal_smoke[present3])

p3d <- plot_ly(
  df3,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~smoke_label, colors = cols3,
  marker = list(size = 6, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text",
  showlegend = TRUE
) %>%
  layout(
    title = "t-SNE 3D - kids (smoking)",
    scene = list(
      xaxis = list(title = "t-SNE 1"),
      yaxis = list(title = "t-SNE 2"),
      zaxis = list(title = "t-SNE 3"),
      aspectmode = "data"
    ),
    legend = list(title = list(text = "Smoking"))
  )

p3d
htmlwidgets::saveWidget(p3d, "results/tsne/tsne3d-09.html",
  selfcontained = TRUE
)


################################################################################
############################# PLOT 10 - parity kids only #######################
################################################################################
# pick the parity column in df
parity_col_df <- if ("parity" %in% names(df)) {
  "parity"
} else if ("mother_birthcard_parity" %in% names(df)) {
  "mother_birthcard_parity"
} else {
  stop("No parity column found in df.")
}

df_par <- df %>%
  mutate(
    parity_num = suppressWarnings(as.numeric(.data[[parity_col_df]])),
    parity_label = case_when(
      is.na(parity_num) ~ "Unknown",
      parity_num == 0 ~ "0",
      parity_num == 1 ~ "1",
      parity_num == 2 ~ "2",
      parity_num >= 3 ~ "3+"
    ),
    parity_label = fct_relevel(parity_label, "0", "1", "2", "3+", "Unknown")
  )

pal_par <- c(
  "0" = "#56B4E9", "1" = "#E69F00", "2" = "#009E73", "3+" = "#CC79A7",
  "Unknown" = "grey70"
)
present <- levels(droplevels(df_par$parity_label))

p2d <- ggplot(df_par, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = parity_label), size = 2, alpha = 0.9) +
  scale_color_manual(
    values = pal_par[present], breaks = present,
    name = "Parity"
  ) +
  labs(
    title = "t-SNE 2D - kids (parity)",
    x = "t-SNE Dimension 1", y = "t-SNE Dimension 2"
  ) +
  theme_phip() +
  theme(legend.position = "top", text = element_text(size = 25))

p2d
ggsave("results/tsne/tsne-10.png",
  plot = p2d, bg = "white", dpi = 300, height = 15, width = 15,
  units = "cm"
)

# --- 3D plotly: kids colored by parity ---

# pick the parity column in 3D df
parity_col_df3 <- if ("parity" %in% names(tsne_df_3d_kid)) {
  "parity"
} else if ("mother_birthcard_parity" %in% names(tsne_df_3d_kid)) {
  "mother_birthcard_parity"
} else {
  stop("No parity column found in tsne_df_3d_kid.")
}

df3 <- tsne_df_3d_kid %>%
  mutate(
    parity_num = suppressWarnings(as.numeric(.data[[parity_col_df3]])),
    parity_label = case_when(
      is.na(parity_num) ~ "Unknown",
      parity_num == 0 ~ "0",
      parity_num == 1 ~ "1",
      parity_num == 2 ~ "2",
      parity_num >= 3 ~ "3+"
    ),
    parity_label = fct_relevel(parity_label, "0", "1", "2", "3+", "Unknown"),
    hover_text = paste0(
      "Subject: ", subject_id, "<br>",
      "Dyade: ", dyade, "<br>",
      "Parity: ", parity_label, "<br>",
      "Timepoint: ", timepoint_factor, "<br>",
      "tSNE1: ", sprintf("%.3f", tSNE1), "<br>",
      "tSNE2: ", sprintf("%.3f", tSNE2), "<br>",
      "tSNE3: ", sprintf("%.3f", tSNE3)
    )
  )

present3 <- levels(droplevels(df3$parity_label))
cols3 <- unname(pal_par[present3])

p3d <- plot_ly(
  df3,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~parity_label, colors = cols3,
  marker = list(size = 6, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text",
  showlegend = TRUE
) %>%
  layout(
    title = "t-SNE 3D - kids (parity)",
    scene = list(
      xaxis = list(title = "t-SNE 1"),
      yaxis = list(title = "t-SNE 2"),
      zaxis = list(title = "t-SNE 3"),
      aspectmode = "data"
    ),
    legend = list(title = list(text = "Parity"))
  )

p3d
htmlwidgets::saveWidget(p3d, "results/tsne/tsne3d-10.html",
  selfcontained = TRUE
)

################################################################################
############################# PLOT 11 - delivery_place kids only ###############
################################################################################
df_place <- df %>%
  mutate(
    place_label = case_when(
      is.na(delivery_place) ~ NA_character_,
      tolower(delivery_place) %in% c("hospital", "clinic") ~ "Hospital",
      tolower(delivery_place) %in% c("home", "at home") ~ "Home",
      TRUE ~ "Other"
    ),
    place_label = fct_explicit_na(place_label, na_level = "Unknown"),
    place_label = fct_relevel(
      place_label, "Hospital", "Home", "Other",
      "Unknown"
    )
  )

# color-blind palette
pal_place <- c(
  "Hospital" = "#0072B2",
  "Home"     = "#E69F00",
  "Other"    = "#CC79A7",
  "Unknown"  = "grey70"
)

present <- levels(droplevels(df_place$place_label))

p2d <- ggplot(df_place, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = place_label), size = 2, alpha = 0.9) +
  scale_color_manual(
    values = pal_place[present],
    breaks = present,
    name = "Delivery place"
  ) +
  labs(
    title = "t-SNE 2D - kids (delivery place)",
    x = "t-SNE Dimension 1", y = "t-SNE Dimension 2"
  ) +
  theme_phip() +
  theme(legend.position = "top", text = element_text(size = 25))

p2d
ggsave("results/tsne/tsne-11.png",
  plot = p2d, bg = "white", dpi = 300, height = 15, width = 15, units = "cm"
)

# --- 3D plotly: kids colored by delivery_place ---
df3 <- tsne_df_3d_kid %>%
  mutate(
    place_label = case_when(
      is.na(delivery_place) ~ NA_character_,
      tolower(delivery_place) %in% c("hospital", "clinic") ~ "Hospital",
      tolower(delivery_place) %in% c("home", "at home") ~ "Home",
      TRUE ~ "Other"
    ),
    place_label = fct_explicit_na(place_label, na_level = "Unknown"),
    place_label = fct_relevel(
      place_label, "Hospital", "Home", "Other",
      "Unknown"
    ),
    hover_text = paste0(
      "Subject: ", subject_id, "<br>",
      "Dyade: ", dyade, "<br>",
      "Delivery place: ", place_label, "<br>",
      "Timepoint: ", timepoint_factor, "<br>",
      "tSNE1: ", sprintf("%.3f", tSNE1), "<br>",
      "tSNE2: ", sprintf("%.3f", tSNE2), "<br>",
      "tSNE3: ", sprintf("%.3f", tSNE3)
    )
  )

present3 <- levels(droplevels(df3$place_label))
cols3 <- unname(pal_place[present3])

p3d <- plot_ly(
  df3,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~place_label, colors = cols3,
  marker = list(size = 6, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text",
  showlegend = TRUE
) %>%
  layout(
    title = "t-SNE 3D - kids (delivery place)",
    scene = list(
      xaxis = list(title = "t-SNE 1"),
      yaxis = list(title = "t-SNE 2"),
      zaxis = list(title = "t-SNE 3"),
      aspectmode = "data"
    ),
    legend = list(title = list(text = "Delivery place"))
  )

p3d
htmlwidgets::saveWidget(p3d, "results/tsne/tsne3d-11.html",
  selfcontained = TRUE
)


################################################################################
############################# ONLY MONTH 12 tSNEs ##############################
################################################################################

## 0) Use a matrix either way
Dm <- if (inherits(D, "dist")) as.matrix(D) else D
stopifnot(identical(rownames(Dm), colnames(Dm))) # must be square

## 1) Parse names and index kid_serum @ T8
nms <- rownames(Dm)
parts <- str_split_fixed(nms, "__", 4)
colnames(parts) <- c("subject_id", "big_group", "timepoint_factor", "dyade")
meta <- as.data.frame(parts, stringsAsFactors = FALSE)

idx_T8 <- meta$big_group == "kid_serum" & meta$timepoint_factor == "T8"

D_kid_T8 <- Dm[idx_T8, idx_T8, drop = FALSE]

parsed_kid_T8 <- tibble(
  id_full          = nms[idx_T8],
  subject_id       = parts[idx_T8, 1],
  big_group        = parts[idx_T8, 2],
  timepoint_factor = parts[idx_T8, 3],
  dyade            = parts[idx_T8, 4]
)

## 2) t-SNE (safe perplexity for subset size)
set.seed(15092025)
n_T8 <- nrow(D_kid_T8)
perp_T8 <- max(5, min(30, floor((n_T8 - 1) / 3))) # Rtsne rule: perp < n/3

ts2_kid_T8 <- Rtsne(
  D_kid_T8,
  dims = 2, is_distance = TRUE,
  perplexity = perp_T8,
  check_duplicates = FALSE, theta = 0.5
)

ts3_kid_T8 <- Rtsne(
  D_kid_T8,
  dims = 3, is_distance = TRUE,
  perplexity = perp_T8,
  check_duplicates = FALSE, theta = 0.5
)

## 3) Bind coords + metadata (avoid dyade clash)
tsne_df_2d_kid_T8 <- parsed_kid_T8 %>%
  mutate(
    tSNE1 = ts2_kid_T8$Y[, 1],
    tSNE2 = ts2_kid_T8$Y[, 2]
  ) %>%
  left_join(select(meta_kids, -dyade), by = "subject_id")

tsne_df_3d_kid_T8 <- parsed_kid_T8 %>%
  mutate(
    tSNE1 = ts3_kid_T8$Y[, 1],
    tSNE2 = ts3_kid_T8$Y[, 2],
    tSNE3 = ts3_kid_T8$Y[, 3]
  ) %>%
  left_join(select(meta_kids, -dyade), by = "subject_id")


# T8 kids with metadata already joined:
df_T8 <- tsne_df_2d_kid_T8
df3_T8 <- tsne_df_3d_kid_T8
## Assumes you've set:
## df_T8  <- tsne_df_2d_kid_T8
## df3_T8 <- tsne_df_3d_kid_T8

################################################################################
############################# PLOT 12 - feedmode_birth (T8) ####################
################################################################################
df_feed_T8 <- df_T8 %>%
  mutate(
    feed_label = recode(feedmode_birth,
      "BF" = "Breastfeeding",
      "FF" = "Formula",
      "FF/BF" = "Mixed",
      .default = NA_character_
    ),
    feed_label = fct_explicit_na(feed_label, na_level = "Unknown"),
    feed_label = fct_relevel(
      feed_label, "Breastfeeding", "Mixed", "Formula",
      "Unknown"
    )
  )

pal_feed <- c(
  "Breastfeeding" = "#009E73", "Mixed" = "#E69F00", "Formula" = "#CC79A7",
  "Unknown" = "grey70"
)
present <- levels(droplevels(df_feed_T8$feed_label))

p12 <- ggplot(df_feed_T8, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = feed_label), size = 2, alpha = 1) +
  scale_color_manual(
    values = pal_feed[present], breaks = present,
    name = "Feeding at birth"
  ) +
  labs(
    title = "t-SNE 2D - kids T8 (feeding at birth)", x = "t-SNE 1",
    y = "t-SNE 2"
  ) +
  theme_phip() +
  theme(legend.position = "top", text = element_text(size = 25))
ggsave("results/tsne/tsne-12.png",
  plot = p12, bg = "white", dpi = 300,
  height = 15, width = 15, units = "cm"
)

df3_feed_T8 <- df3_T8 %>%
  mutate(
    feed_label = recode(feedmode_birth,
      "BF" = "Breastfeeding", "FF" = "Formula", "FF/BF" = "Mixed",
      .default = NA_character_
    ),
    feed_label = fct_explicit_na(feed_label, na_level = "Unknown"),
    feed_label = fct_relevel(
      feed_label, "Breastfeeding", "Mixed", "Formula",
      "Unknown"
    ),
    hover_text = paste0(
      "Subject: ", subject_id, "<br>",
      "Dyade: ", dyade, "<br>",
      "Feeding at birth: ", feed_label, "<br>",
      "tSNE1: ", sprintf("%.3f", tSNE1), "<br>",
      "tSNE2: ", sprintf("%.3f", tSNE2), "<br>",
      "tSNE3: ", sprintf("%.3f", tSNE3)
    )
  )
present3 <- levels(droplevels(df3_feed_T8$feed_label))
cols3 <- unname(pal_feed[present3])

p12_3d <- plot_ly(
  df3_feed_T8,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~feed_label, colors = cols3,
  marker = list(size = 6, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text", showlegend = TRUE
) %>% layout(
  title = "t-SNE 3D - kids T8 (feeding at birth)",
  scene = list(
    xaxis = list(title = "t-SNE 1"), yaxis = list(title = "t-SNE 2"),
    zaxis = list(title = "t-SNE 3"), aspectmode = "data"
  ),
  legend = list(title = list(text = "Feeding at birth"))
)
saveWidget(p12_3d, "results/tsne/tsne3d-12.html", selfcontained = TRUE)


################################################################################
############################# PLOT 13 - infant sex (T8) ########################
################################################################################
df_sex_T8 <- df_T8 %>%
  mutate(
    sex_label = recode(infant_sex,
      male = "Male", female = "Female",
      .default = NA_character_
    ),
    sex_label = fct_explicit_na(sex_label, na_level = "Unknown"),
    sex_label = fct_relevel(sex_label, "Female", "Male", "Unknown")
  )

pal_sex <- c("Female" = "#CC79A7", "Male" = "#56B4E9", "Unknown" = "grey70")
present <- levels(droplevels(df_sex_T8$sex_label))

p13 <- ggplot(df_sex_T8, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = sex_label), size = 2, alpha = 0.9) +
  scale_color_manual(
    values = pal_sex[present], breaks = present,
    name = "Infant sex"
  ) +
  labs(
    title = "t-SNE 2D - kids T8 (by infant sex)", x = "t-SNE 1",
    y = "t-SNE 2"
  ) +
  theme_phip() +
  theme(legend.position = "top", text = element_text(size = 25))
ggsave("results/tsne/tsne-13.png",
  plot = p13, bg = "white", dpi = 300,
  height = 15, width = 15, units = "cm"
)

df3_sex_T8 <- df3_T8 %>%
  mutate(
    sex_label = recode(infant_sex,
      male = "Male", female = "Female",
      .default = NA_character_
    ),
    sex_label = fct_explicit_na(sex_label, na_level = "Unknown"),
    sex_label = fct_relevel(sex_label, "Female", "Male", "Unknown"),
    hover_text = paste0(
      "Subject: ", subject_id, "<br>",
      "Dyade: ", dyade, "<br>",
      "Infant sex: ", sex_label, "<br>",
      "tSNE1: ", sprintf("%.3f", tSNE1), "<br>",
      "tSNE2: ", sprintf("%.3f", tSNE2), "<br>",
      "tSNE3: ", sprintf("%.3f", tSNE3)
    )
  )
present3 <- levels(droplevels(df3_sex_T8$sex_label))
cols3 <- unname(pal_sex[present3])

p13_3d <- plot_ly(
  df3_sex_T8,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~sex_label, colors = cols3,
  marker = list(size = 6, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text", showlegend = TRUE
) %>% layout(
  title = "t-SNE 3D - kids T8 (by infant sex)",
  scene = list(
    xaxis = list(title = "t-SNE 1"), yaxis = list(title = "t-SNE 2"),
    zaxis = list(title = "t-SNE 3"), aspectmode = "data"
  ),
  legend = list(title = list(text = "Infant sex"))
)
saveWidget(p13_3d, "results/tsne/tsne3d-13.html", selfcontained = TRUE)


################################################################################
############################# PLOT 14 - feedmode at 3m (T8) ####################
################################################################################
df_m3_T8 <- df_T8 %>%
  mutate(
    m3_label = case_when(
      feedmode_m3 == "BF" ~ "Breastfeeding",
      feedmode_m3 == "FF" ~ "Formula",
      feedmode_m3 == "FF/BF" ~ "Mixed",
      grepl("^Otherwise", feedmode_m3 %||% "") ~ "Other",
      TRUE ~ NA_character_
    ),
    m3_label = fct_explicit_na(m3_label, na_level = "Unknown"),
    m3_label = fct_relevel(
      m3_label, "Breastfeeding", "Mixed", "Formula",
      "Other", "Unknown"
    )
  )

pal_m3 <- c(
  "Breastfeeding" = "#009E73", "Mixed" = "#E69F00",
  "Formula" = "#CC79A7", "Other" = "#9467BD", "Unknown" = "grey70"
)
present <- levels(droplevels(df_m3_T8$m3_label))

p14 <- ggplot(df_m3_T8, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = m3_label), size = 2, alpha = 0.9) +
  scale_color_manual(
    values = pal_m3[present], breaks = present,
    name = "Feeding at 3 mo"
  ) +
  labs(
    title = "t-SNE 2D - kids T8 (feeding at 3 months)",
    x = "t-SNE 1", y = "t-SNE 2"
  ) +
  theme_phip() +
  theme(legend.position = "top", text = element_text(size = 25))
ggsave("results/tsne/tsne-14.png",
  plot = p14, bg = "white", dpi = 300,
  height = 15, width = 15, units = "cm"
)

df3_m3_T8 <- df3_T8 %>%
  mutate(
    m3_label = case_when(
      feedmode_m3 == "BF" ~ "Breastfeeding",
      feedmode_m3 == "FF" ~ "Formula",
      feedmode_m3 == "FF/BF" ~ "Mixed",
      grepl("^Otherwise", feedmode_m3 %||% "") ~ "Other",
      TRUE ~ NA_character_
    ),
    m3_label = fct_explicit_na(m3_label, na_level = "Unknown"),
    m3_label = fct_relevel(
      m3_label, "Breastfeeding", "Mixed", "Formula",
      "Other", "Unknown"
    ),
    hover_text = paste0(
      "Subject: ", subject_id, "<br>",
      "Dyade: ", dyade, "<br>",
      "Feeding @3m: ", m3_label, "<br>",
      "tSNE1: ", sprintf("%.3f", tSNE1), "<br>",
      "tSNE2: ", sprintf("%.3f", tSNE2), "<br>",
      "tSNE3: ", sprintf("%.3f", tSNE3)
    )
  )
present3 <- levels(droplevels(df3_m3_T8$m3_label))
cols3 <- unname(pal_m3[present3])

p14_3d <- plot_ly(
  df3_m3_T8,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~m3_label, colors = cols3,
  marker = list(size = 6, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text", showlegend = TRUE
) %>% layout(
  title = "t-SNE 3D - kids T8 (feeding at 3 months)",
  scene = list(
    xaxis = list(title = "t-SNE 1"), yaxis = list(title = "t-SNE 2"),
    zaxis = list(title = "t-SNE 3"), aspectmode = "data"
  ),
  legend = list(title = list(text = "Feeding at 3 mo"))
)
saveWidget(p14_3d, "results/tsne/tsne3d-14.html", selfcontained = TRUE)


################################################################################
############################# PLOT 15 - smoking (T8) ###########################
################################################################################
df_smoke_T8 <- df_T8 %>%
  mutate(
    smoking = smoking, # keep if already present
    smoke_label = case_when(
      is.na(smoking) ~ NA_character_,
      tolower(smoking) %in% c("yes", "y", "1", "true") ~ "Yes",
      tolower(smoking) %in% c("no", "n", "0", "false") ~ "No",
      TRUE ~ "Other"
    ),
    smoke_label = fct_explicit_na(smoke_label, na_level = "Unknown"),
    smoke_label = fct_relevel(smoke_label, "No", "Yes", "Unknown", "Other")
  )

pal_smoke <- c(
  "No" = "#009E73", "Yes" = "#D55E00", "Unknown" = "grey70",
  "Other" = "#9467BD"
)
present <- levels(droplevels(df_smoke_T8$smoke_label))

p15 <- ggplot(df_smoke_T8, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = smoke_label), size = 2, alpha = 0.9) +
  scale_color_manual(
    values = pal_smoke[present], breaks = present,
    name = "Smoking"
  ) +
  labs(title = "t-SNE 2D - kids T8 (smoking)", x = "t-SNE 1", y = "t-SNE 2") +
  theme_phip() +
  theme(legend.position = "top", text = element_text(size = 25))
ggsave("results/tsne/tsne-15.png",
  plot = p15, bg = "white", dpi = 300,
  height = 15, width = 15, units = "cm"
)

df3_smoke_T8 <- df3_T8 %>%
  mutate(
    smoke_label = case_when(
      is.na(smoking) ~ NA_character_,
      tolower(smoking) %in% c("yes", "y", "1", "true") ~ "Yes",
      tolower(smoking) %in% c("no", "n", "0", "false") ~ "No",
      TRUE ~ "Other"
    ),
    smoke_label = fct_explicit_na(smoke_label, na_level = "Unknown"),
    smoke_label = fct_relevel(smoke_label, "No", "Yes", "Unknown", "Other"),
    hover_text = paste0(
      "Subject: ", subject_id, "<br>",
      "Dyade: ", dyade, "<br>",
      "Smoking: ", smoke_label, "<br>",
      "tSNE1: ", sprintf("%.3f", tSNE1), "<br>",
      "tSNE2: ", sprintf("%.3f", tSNE2), "<br>",
      "tSNE3: ", sprintf("%.3f", tSNE3)
    )
  )
present3 <- levels(droplevels(df3_smoke_T8$smoke_label))
cols3 <- unname(pal_smoke[present3])

p15_3d <- plot_ly(
  df3_smoke_T8,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~smoke_label, colors = cols3,
  marker = list(size = 6, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text", showlegend = TRUE
) %>% layout(
  title = "t-SNE 3D - kids T8 (smoking)",
  scene = list(
    xaxis = list(title = "t-SNE 1"), yaxis = list(title = "t-SNE 2"),
    zaxis = list(title = "t-SNE 3"), aspectmode = "data"
  ),
  legend = list(title = list(text = "Smoking"))
)
saveWidget(p15_3d, "results/tsne/tsne3d-15.html", selfcontained = TRUE)

################################################################################
############################# PLOT 16 - parity (T8) ############################
################################################################################
parity_col_df_T8 <- if ("parity" %in% names(df_T8)) {
  "parity"
} else if ("mother_birthcard_parity" %in% names(df_T8)) {
  "mother_birthcard_parity"
} else {
  stop("No parity column found in df_T8.")
}

df_par_T8 <- df_T8 %>%
  mutate(
    parity_num = suppressWarnings(as.numeric(.data[[parity_col_df_T8]])),
    parity_label = case_when(
      is.na(parity_num) ~ "Unknown",
      parity_num == 0 ~ "0",
      parity_num == 1 ~ "1",
      parity_num == 2 ~ "2",
      parity_num >= 3 ~ "3+"
    ),
    parity_label = fct_relevel(parity_label, "0", "1", "2", "3+", "Unknown")
  )

pal_par <- c(
  "0" = "#56B4E9", "1" = "#E69F00", "2" = "#009E73", "3+" = "#CC79A7",
  "Unknown" = "grey70"
)
present <- levels(droplevels(df_par_T8$parity_label))

p16 <- ggplot(df_par_T8, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = parity_label), size = 2, alpha = 0.9) +
  scale_color_manual(
    values = pal_par[present], breaks = present,
    name = "Parity"
  ) +
  labs(title = "t-SNE 2D - kids T8 (parity)", x = "t-SNE 1", y = "t-SNE 2") +
  theme_phip() +
  theme(legend.position = "top", text = element_text(size = 25))
ggsave("results/tsne/tsne-16.png",
  plot = p16, bg = "white", dpi = 300,
  height = 15, width = 15, units = "cm"
)

parity_col_df3_T8 <- if ("parity" %in% names(df3_T8)) {
  "parity"
} else if ("mother_birthcard_parity" %in% names(df3_T8)) {
  "mother_birthcard_parity"
} else {
  stop("No parity column found in df3_T8.")
}

df3_par_T8 <- df3_T8 %>%
  mutate(
    parity_num = suppressWarnings(as.numeric(.data[[parity_col_df3_T8]])),
    parity_label = case_when(
      is.na(parity_num) ~ "Unknown",
      parity_num == 0 ~ "0",
      parity_num == 1 ~ "1",
      parity_num == 2 ~ "2",
      parity_num >= 3 ~ "3+"
    ),
    parity_label = fct_relevel(parity_label, "0", "1", "2", "3+", "Unknown"),
    hover_text = paste0(
      "Subject: ", subject_id, "<br>",
      "Dyade: ", dyade, "<br>",
      "Parity: ", parity_label, "<br>",
      "tSNE1: ", sprintf("%.3f", tSNE1), "<br>",
      "tSNE2: ", sprintf("%.3f", tSNE2), "<br>",
      "tSNE3: ", sprintf("%.3f", tSNE3)
    )
  )
present3 <- levels(droplevels(df3_par_T8$parity_label))
cols3 <- unname(pal_par[present3])

p16_3d <- plot_ly(
  df3_par_T8,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~parity_label, colors = cols3,
  marker = list(size = 6, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text", showlegend = TRUE
) %>% layout(
  title = "t-SNE 3D - kids T8 (parity)",
  scene = list(
    xaxis = list(title = "t-SNE 1"), yaxis = list(title = "t-SNE 2"),
    zaxis = list(title = "t-SNE 3"), aspectmode = "data"
  ),
  legend = list(title = list(text = "Parity"))
)
saveWidget(p16_3d, "results/tsne/tsne3d-16.html", selfcontained = TRUE)


################################################################################
############################# PLOT 17 - delivery place (T8) ####################
################################################################################
df_place_T8 <- df_T8 %>%
  mutate(
    place_label = case_when(
      is.na(delivery_place) ~ NA_character_,
      tolower(delivery_place) %in% c("hospital", "clinic") ~ "Hospital",
      tolower(delivery_place) %in% c("home", "at home") ~ "Home",
      TRUE ~ "Other"
    ),
    place_label = fct_explicit_na(place_label, na_level = "Unknown"),
    place_label = fct_relevel(
      place_label, "Hospital", "Home", "Other",
      "Unknown"
    )
  )

pal_place <- c(
  "Hospital" = "#0072B2", "Home" = "#E69F00", "Other" = "#CC79A7",
  "Unknown" = "grey70"
)
present <- levels(droplevels(df_place_T8$place_label))

p17 <- ggplot(df_place_T8, aes(tSNE1, tSNE2)) +
  geom_point(aes(color = place_label), size = 2, alpha = 0.9) +
  scale_color_manual(
    values = pal_place[present], breaks = present,
    name = "Delivery place"
  ) +
  labs(
    title = "t-SNE 2D - kids T8 (delivery place)", x = "t-SNE 1",
    y = "t-SNE 2"
  ) +
  theme_phip() +
  theme(legend.position = "top", text = element_text(size = 25))
ggsave("results/tsne/tsne-17.png",
  plot = p17, bg = "white", dpi = 300,
  height = 15, width = 15, units = "cm"
)

df3_place_T8 <- df3_T8 %>%
  mutate(
    place_label = case_when(
      is.na(delivery_place) ~ NA_character_,
      tolower(delivery_place) %in% c("hospital", "clinic") ~ "Hospital",
      tolower(delivery_place) %in% c("home", "at home") ~ "Home",
      TRUE ~ "Other"
    ),
    place_label = fct_explicit_na(place_label, na_level = "Unknown"),
    place_label = fct_relevel(
      place_label, "Hospital", "Home", "Other",
      "Unknown"
    ),
    hover_text = paste0(
      "Subject: ", subject_id, "<br>",
      "Dyade: ", dyade, "<br>",
      "Delivery place: ", place_label, "<br>",
      "tSNE1: ", sprintf("%.3f", tSNE1), "<br>",
      "tSNE2: ", sprintf("%.3f", tSNE2), "<br>",
      "tSNE3: ", sprintf("%.3f", tSNE3)
    )
  )
present3 <- levels(droplevels(df3_place_T8$place_label))
cols3 <- unname(pal_place[present3])

p17_3d <- plot_ly(
  df3_place_T8,
  x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
  type = "scatter3d", mode = "markers",
  color = ~place_label, colors = cols3,
  marker = list(size = 6, opacity = 0.9),
  text = ~hover_text, hoverinfo = "text", showlegend = TRUE
) %>% layout(
  title = "t-SNE 3D - kids T8 (delivery place)",
  scene = list(
    xaxis = list(title = "t-SNE 1"), yaxis = list(title = "t-SNE 2"),
    zaxis = list(title = "t-SNE 3"), aspectmode = "data"
  ),
  legend = list(title = list(text = "Delivery place"))
)
saveWidget(p17_3d, "results/tsne/tsne3d-17.html", selfcontained = TRUE)

## remove vars and clean
rm(
  list = setdiff(
    ls(envir = .GlobalEnv, all.names = TRUE),
    c("ps_merged_box_bin", "seed_before")
  ),
  envir = .GlobalEnv
)

gc()
gc()
