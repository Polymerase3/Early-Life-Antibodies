################################################################################
########################### COMPUTING STABILITIES ##############################
################################################################################
# this is the main workhorse for the stabilities generally - it computes
# probably everything one could think of and stores it in a single table; so the
# stabilities/correlations between the subjects at a given timpeoint and also
# between the timepoints and correlations between the different groups (mom_serum,
# kid_serum, mom_milk) too. --> the only tricky part is to extract the data from
# the table, as you have to filter it down to the groups/timepoints of interest
con <- ps_merged_box_bin$meta$con

# 1) compute stability using dyade_recoded and keep only kid_serum vs kid_serum
# rows
stab_recoded_tbl <- compute_repertoire_similarity(
  ps_merged_box_bin,
  group_col = "big_group",
  time_col = "timepoint_recoded",
  similarity = "kulczynski",
  mode = "all",
  dyade_col = "dyade_recoded",
  time_mode = "pairwise",
  time_pairing = c("same", "cross"),
  subject_mode = "all",
  drop_self_time = FALSE,
  overwrite = TRUE,
  verbose = FALSE
)

# materialise kid vs kid subset as persistent duckdb table "stab_kid_kid"
if (DBI::dbExistsTable(con, "stab_kid_kid")) {
  DBI::dbExecute(con, "DROP TABLE IF EXISTS stab_kid_kid")
}

# select only kid vs kid rows
stab_recoded_tbl %>%
  filter(group_left == "kid_serum", group_right == "kid_serum") %>%
  compute(name = "stab_kid_kid", temporary = FALSE)

# 2) compute full stability using original dyade (not dyade_recoded)
stab_orig_tbl <- compute_repertoire_similarity(
  ps_merged_box_bin,
  group_col = "big_group",
  time_col = "timepoint_recoded",
  similarity = "kulczynski",
  mode = "all",
  dyade_col = "dyade",
  time_mode = "pairwise",
  time_pairing = c("same", "cross"),
  subject_mode = "all",
  drop_self_time = FALSE,
  overwrite = TRUE,
  verbose = FALSE
)

# 3) remove kid_serum vs kid_serum rows from the dyade-based full table and
# union with the saved kid-kid subset; materialise the final result as
# persistent duckdb table named "stab_df_full"
if (DBI::dbExistsTable(con, "stab_df_full")) {
  DBI::dbExecute(con, "DROP TABLE IF EXISTS stab_df_full")
}

# drop kid-kid from dyade-based result and append kid-kid from recoded result
stab_orig_tbl %>%
  filter(!(group_left == "kid_serum" & group_right == "kid_serum")) %>%
  union_all(tbl(con, "stab_kid_kid")) %>%
  compute(name = "stab_df_full", temporary = FALSE)

# 4) expose final duckdb table as an R tbl named stab_df_full
stab_df_full <- dplyr::tbl(con, "stab_df_full")

# 5) recoding the time_left and time_right back to the usual, descriptive labels
# small mapping tibble
con <- dbplyr::remote_con(stab_df_full)
tbl_name <- dbplyr::remote_name(stab_df_full)

# map T --> P/M
map_t2p <- tibble::tibble(
  from = c("T0", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8"),
  to   = c("P12", "P28", "M0", "M1", "M1", "M3", "M3", "M6", "M12")
)

# temp table with map
copy_to(con, map_t2p, name = "map_t2p_tmp", temporary = TRUE, overwrite = TRUE)

# update the solumns in stab_df_full
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

# ------------------ recoding the dyades ---------------------------------------
# thomas wanted the siblings and twins "more in the middle" and not at the
# borders of the plot. It was tricky, cause in my phip_data object, one twins
# had probably a large dyade variable value and as i sort the observations in
# the heatmaps based on dyades, they always landed on the right end of the plot.
# I figured out the solution would actually be to swap the dyade numbers with
# some other dyades that were more in the middle --> it doesnt really matter
# where they are, as the order of the samples in the heatmap is totally
# arbitrary

# work on a copy for now to be safe
stab_df_swapped <- stab_df_full

# 1) get the max dyad without pulling the whole table
max_dyad <- stab_df_full %>%
  transmute(
    dyad_left  = as.integer(dyad_left),
    dyad_right = as.integer(dyad_right)
  ) %>%
  summarise(
    m1 = max(dyad_left, na.rm = TRUE),
    m2 = max(dyad_right, na.rm = TRUE)
  ) %>%
  collect() %>%
  {
    max(c(.$m1, .$m2), na.rm = TRUE)
  }

# 2) define the bands (bijective, 40 ids each)
hi_range <- seq(max_dyad - 39L, max_dyad) # the last 40 dyads
mid_range <- 40L:79L # the 40–79 band

# 3) build the mapping (as character to match your columns)
map_df <- tibble::tibble(
  old = as.character(c(hi_range, mid_range)),
  new = as.character(c(mid_range, hi_range))
)

# 4) apply the mapping to both dyad columns via left joins (lazy-safe)
stab_df_swapped <- stab_df_full %>%
  # swap dyad_left
  left_join(map_df, by = c("dyad_left" = "old"), copy = TRUE) %>%
  mutate(dyad_left = dplyr::coalesce(new, dyad_left)) %>%
  select(-new) %>%
  # swap dyad_right
  left_join(map_df, by = c("dyad_right" = "old"), copy = TRUE) %>%
  mutate(dyad_right = dplyr::coalesce(new, dyad_right)) %>%
  select(-new)

# overwrite the original table
stab_df_full <- stab_df_swapped

################################################################################
############################## PLOTTING ########################################
################################################################################
# ----------------------------- helpers ----------------------------------------
# reader-friendly labels
pretty_group <- function(g) gsub("_", " ", g) |> tools::toTitleCase()
pretty_time <- function(t) {
  if (t == "M0") {
    return("0 months")
  }
  if (grepl("^M\\d+$", t)) {
    n <- as.integer(sub("^M", "", t))
    if (n == 1) "1 month" else sprintf("%d months", n)
  } else if (grepl("^P\\d+$", t)) {
    n <- as.integer(sub("^P", "", t))
    if (n == 1) "1 week (pregnancy)" else sprintf("%d weeks (pregnancy)", n)
  } else {
    t
  }
}

# safe slug for filenames
slug <- function(...) {
  x <- paste(..., sep = "__") # single, internal sep
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  tolower(trimws(x, which = "both", whitespace = "_"))
}

# ------------------------ enumerate all pairs ---------------------------------
# use only states that actually exist in the table (left and right sides
# separately)
left_states <- stab_df_full %>%
  distinct(g1 = group_left, t1 = time_left) %>%
  collect()
right_states <- stab_df_full %>%
  distinct(g2 = group_right, t2 = time_right) %>%
  collect()

# cartesian product of left * right (all possible pairs present in data)
pairs_df <- tidyr::crossing(left_states, right_states)

# deterministic ordering by group, then by time
time_order <- c("P12", "P28", "M0", "M1", "M2", "M3", "M6", "M12")
ord_time <- function(x) {
  factor(x,
    levels = unique(c(time_order, setdiff(x, time_order))),
    ordered = TRUE
  )
}

pairs_df <- pairs_df %>%
  mutate(
    g1o = as.character(g1), g2o = as.character(g2),
    t1o = ord_time(t1), t2o = ord_time(t2)
  ) %>%
  arrange(g1o, t1o, g2o, t2o) %>%
  select(-g1o, -g2o, -t1o, -t2o)

# ----------------------- output dir & numbering -------------------------------
dir.create("results/heatmaps_tables", recursive = TRUE, showWarnings = FALSE)
n <- nrow(pairs_df)
pad <- max(2, nchar(as.character(n))) # zero-padding (min 2 digits)

# log results
export_log <- vector("list", n)
cat(sprintf("Generating %d heatmap(s)...\n", n))

# --- helper: which subject IDs to drop for a given comparison ---
excluded_subjects <- function(g1, t1, g2, t2) {
  ex <- character(0)

  # wait for the decision about outliers!

  # specific heatmaps (both orders)
  # if ((g1 == "mom_serum" && t1 == "M0" && g2 == "mom_milk" && t2 == "M1") ||
  #     (g2 == "mom_serum" && t2 == "M0" && g1 == "mom_milk" && t1 == "M1")) {
  #   ex <- c(ex, "011783")
  # }
  #
  # if ((g1 == "mom_serum" && t1 == "M0" && g2 == "mom_milk" && t2 == "M12") ||
  #     (g2 == "mom_serum" && t2 == "M0" && g1 == "mom_milk" && t1 == "M12")) {
  #   ex <- c(ex, "008048")
  # }
  #
  # # group-level rules (any timepoints, both orders)
  # if ((g1 == "mom_milk" && g2 == "mom_milk")) {
  #   ex <- c(ex, "200736", "011783")
  # }
  #
  # if ((g1 == "mom_milk" && g2 == "kid_serum") ||
  #     (g2 == "mom_milk" && g1 == "kid_serum")) {
  #   ex <- c(ex, "003520")
  # }

  unique(ex)
}

# ---------------------------- main loop ---------------------------------------
for (i in seq_len(n)) {
  g1 <- pairs_df$g1[[i]]
  t1 <- pairs_df$t1[[i]]
  g2 <- pairs_df$g2[[i]]
  t2 <- pairs_df$t2[[i]]

  # decide pairing/dyads and extra args
  pairing <- if (identical(t1, t2)) "same" else "cross"

  dyads <- if (pairing == "same") {
    # same timepoint: keep dyad order, preserve diagonal
    "order"
  } else if (t1 == t2) {
    # cross time within the same group: do not enforce dyads
    "no"
  } else {
    # cross time & cross group: keep ordered dyads
    "order"
  }

  # extra args only for cross comparisons
  extra <- if (pairing == "cross") {
    list(grid = "intersection", dedupe = "mean")
  } else {
    list()
  }

  # file name: heat-XX-<g1>-<t1>__<g2>-<t2>.svg
  file_stub <- sprintf(paste0("heat-%0", pad, "d-"), i)
  file_name <- paste0(
    file_stub,
    slug(paste0(g1, "-", t1), paste0(g2, "-", t2)),
    "_inferno_final", # palette name
    ".svg"
  )

  outfile <- file.path("results/heatmaps_tables", file_name)

  status <- "ok"
  err <- NA_character_

  # label helpers
  wrap_title <- function(s, max = 45) {
    if (is.null(s) || is.na(s)) {
      return(s)
    }
    if (nchar(s) > max && grepl(" vs ", s, fixed = TRUE)) {
      return(sub(" vs ", "\nvs ", s, fixed = TRUE))
    }
    str_wrap(s, width = max)
  }

  wrap_lab <- function(s, width = 45) {
    if (is.null(s) || is.na(s)) {
      return(s)
    }
    str_wrap(s, width = width)
  }

  xlab <- sprintf("%s - %s", pretty_time(t2), pretty_group(g2))
  ylab <- sprintf("%s - %s", pretty_time(t1), pretty_group(g1))

  title_raw <- if (pairing == "same") {
    sprintf(
      "%s - %s vs %s", pretty_time(t1), pretty_group(g1),
      pretty_group(g2)
    )
  } else if (g1 == g2) {
    sprintf(
      "%s vs %s - %s", pretty_time(t1), pretty_time(t2),
      pretty_group(g1)
    )
  } else {
    sprintf(
      "%s - %s vs %s - %s",
      pretty_time(t1), pretty_group(g1),
      pretty_time(t2), pretty_group(g2)
    )
  }

  pal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

  # build a filtered (still-lazy) sim table for this comparison; doesnt have to
  # pull/collect everything, which is leaner
  ids_to_drop <- excluded_subjects(g1, t1, g2, t2)

  sim_tbl_i <- if (length(ids_to_drop)) {
    stab_df_full %>%
      dplyr::filter(!(subject_left %in% ids_to_drop |
        subject_right %in% ids_to_drop))
  } else {
    stab_df_full
  }

  # wrap plotting in tryCatch so a failure doesnt stop the loop
  p <- tryCatch(
    {
      do.call(
        plot_similarity_heatmap,
        c(list(
          sim_tbl = stab_df_full,
          groups  = c(g1, g2),
          times   = c(t1, t2),
          pairing = pairing,
          dyads   = dyads,
          cluster = "none"
        ), extra)
      ) +
        # scale_fill_gradientn(
        #   colours = pal,
        #   values  = c(0, 0.25, 0.5, 0.75, 1),  # keeps yellow centered at 0.5
        #   limits  = c(0,1),
        #   oob     = scales::squish,
        #   name    = "Kulczynski similarity"
        # ) +
        scale_fill_viridis_c(
          option = "inferno",
          direction = -1, # use -1 to invert
          limits = c(0, 1),
          oob = scales::squish,
          name = "Kulczynski similarity"
        ) +
        labs(
          x = wrap_lab(xlab, 40),
          y = wrap_lab(ylab, 40),
          title = wrap_title(title_raw, 42)
        ) +
        guides(
          fill = guide_colorbar(
            title.position = "top",
            title.hjust    = 0.5,
            barwidth       = grid::unit(8, "cm"),
            barheight      = grid::unit(0.45, "cm")
          )
        ) +
        # reasonable tick thinning for dense axes
        scale_x_discrete(
          breaks = function(x) {
            idx <- unique(c(1, seq(5, length(x), by = 5)))
            x[idx]
          },
          labels = function(b) {
            if (length(b) <= 1) {
              return(b)
            }
            c(1, 5 * seq_len(length(b) - 1))
          }
        ) +
        scale_y_discrete(
          breaks = function(y) {
            idx <- unique(c(1, seq(5, length(y), by = 5)))
            y[idx]
          },
          labels = function(b) {
            if (length(b) <= 1) {
              return(b)
            }
            c(1, 5 * seq_len(length(b) - 1))
          }
        ) +
        theme(
          legend.position      = "top",
          legend.direction     = "horizontal",
          legend.justification = "center",
          plot.title           = element_text(lineheight = 0.8, size = 16),
          legend.box.margin    = margin(b = 6),
          text                 = element_text(size = 16),
          axis.text            = element_text(size = 12)
        )
    },
    error = function(e) {
      status <<- "plot_error"
      err <<- conditionMessage(e)
      NULL
    }
  )

  # if plotting succeeded, try saving; otherwise skip file write
  if (!is.null(p)) {
    # --- build CSV in the exact on-plot order (visual match) ------------------
    hm_df <- attr(p, "heatmap_df")
    axes <- attr(p, "axis_levels")
    if (!is.null(hm_df) && !is.null(axes)) {
      # pivot to wide and enforce axis order
      wide <- tidyr::pivot_wider(
        hm_df,
        id_cols     = subject_left,
        names_from  = subject_right,
        values_from = similarity
      )

      row_levels_plot <- axes$rows # bottom -> top in ggplot
      col_levels_plot <- axes$cols # left   -> right in ggplot

      # Make top row in CSV = top of the plot  (reverse the ggplot row order)
      wide <- wide %>%
        dplyr::mutate(subject_left = as.character(subject_left)) %>%
        dplyr::arrange(factor(subject_left, levels = rev(row_levels_plot)))

      # Reorder columns to match X axis left -> right
      keep_cols <- intersect(
        col_levels_plot,
        setdiff(names(wide), "subject_left")
      )
      wide <- wide[, c("subject_left", keep_cols), drop = FALSE]

      # filenames
      outfile_csv <- sub("\\.svg$", ".csv", outfile)
      readr::write_csv(wide, outfile_csv, na = "")
    }

    tryCatch(
      {
        ggsave(outfile,
          plot = p, bg = "white", dpi = 300, height = 15,
          width = 15, units = "cm"
        )
      },
      error = function(e) {
        status <<- "save_error"
        err <<- conditionMessage(e)
      }
    )
  }

  # record log row
  export_log[[i]] <- data.frame(
    idx = i,
    g1 = g1, t1 = t1,
    g2 = g2, t2 = t2,
    pairing = pairing,
    dyads = dyads,
    outfile = outfile,
    status = status,
    error = err,
    stringsAsFactors = FALSE
  )

  # console feedback (short)
  if (identical(status, "ok")) {
    cat(sprintf("[%0*d] OK  -> %s\n", pad, i, basename(outfile)))
  } else {
    cat(sprintf(
      "[%0*d] SKIP (%s): %s\n", pad, i, status,
      if (is.na(err)) "" else err
    ), sep = "")
  }
}

# write a CSV log with all attempts (ok + failures)
log_df <- dplyr::bind_rows(export_log)
readr::write_csv(log_df, file.path(
  "results/heatmaps_tables",
  "heatmap_export_log.csv"
))

cat("Done.\n")

################################################################################
######################### STABILITY BOXPLOT ####################################
################################################################################
# ---------- helper: map printed time labels (m0 -> b) ----------
# function to convert raw time tokens for printing/labels: m0 -> b;
# everything else unchanged
map_time_label <- function(t) {
  if (is.null(t)) {
    return(t)
  }
  dplyr::recode(as.character(t), M0 = "B", .default = as.character(t))
}

# ---------- 0) compute pairwise similarities (phip) ----------
stab_df_full <- compute_repertoire_similarity(
  ps_merged_box_bin,
  group_col      = "big_group",
  time_col       = "timepoint_recoded",
  similarity     = "kulczynski",
  mode           = "all",
  dyade_col      = "dyade_recoded",
  time_mode      = "pairwise",
  time_pairing   = c("same", "cross"),
  subject_mode   = "all",
  drop_self_time = FALSE,
  overwrite      = TRUE,
  verbose        = TRUE
)

# ---------- 1) recode t0..t8 -> p/m labels directly in duckdb ----------
# small mapping table and update left/right time columns in the db table
con <- dbplyr::remote_con(stab_df_full)
tbl_name <- dbplyr::remote_name(stab_df_full)

map_t2p <- tibble::tibble(
  from = c("T0", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8"),
  to   = c("P12", "P28", "M0", "M1", "M1", "M3", "M3", "M6", "M12")
)
copy_to(con, map_t2p, name = "map_t2p_tmp", temporary = TRUE, overwrite = TRUE)

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

# ---------- 2) specification of comparisons (labels/order/blocks) ----------
# note: labels here already use the printed form "b" for birth (m_b / i_b)
cmp_specs <- tribble(
  ~label, ~group_a, ~time_a, ~group_b, ~time_b, ~match_by, ~block,
  # pregnancy (mom_serum)
  "M_P12 vs M_B", "mom_serum", "P12", "mom_serum", "M0", "subject", "Pregnancy",
  "M_P28 vs M_B", "mom_serum", "P28", "mom_serum", "M0", "subject", "Pregnancy",
  "M_P28 vs M_P12", "mom_serum", "P28", "mom_serum", "P12", "subject", "Pregnancy",
  # bridge at birth: mom at birth vs infant at birth (dyad)
  "M_B vs I_B", "mom_serum", "M0", "kid_serum", "M0", "dyad", "Infants",
  # infants (kid_serum)
  "I_B vs I_M3", "kid_serum", "M0", "kid_serum", "M3", "subject", "Infants",
  "I_M3 vs I_M12", "kid_serum", "M3", "kid_serum", "M12", "subject", "Infants",
  "I_B vs I_M12", "kid_serum", "M0", "kid_serum", "M12", "subject", "Infants",
  # breast milk (mom_milk)
  "BM_M1 vs BM_M3", "mom_milk", "M1", "mom_milk", "M3", "subject", "Breast milk",
  "BM_M3 vs BM_M6", "mom_milk", "M3", "mom_milk", "M6", "subject", "Breast milk",
  "BM_M1 vs BM_M6", "mom_milk", "M1", "mom_milk", "M6", "subject", "Breast milk"
)
label_order <- cmp_specs$label

# ---------- 3) normalise left/right; db-safe extractor ----------
norm_pairs <- stab_df_full %>%
  transmute(
    similarity,
    grpA = group_left, tmA = time_left,
    subjA = subject_left, dyadA = dyad_left,
    grpB = group_right, tmB = time_right,
    subjB = subject_right, dyadB = dyad_right
  )

# empty template returned when no matches
tmpl <- tibble(
  label = character(), block = character(), id = character(),
  similarity = double()
)

# helper: safely pull comparison pairs from the db (left/right normalised)
pull_cmp_safe <- function(label, group_a, time_a, group_b, time_b, match_by,
                          block, dat = norm_pairs) {
  gA <- group_a
  tA <- time_a
  gB <- group_b
  tB <- time_b
  by <- match_by

  # 0) gather ids to exclude for this comparison
  ex <- excluded_subjects(gA, tA, gB, tB)

  # 1) match pairs irrespective of side + filter excluded subjects
  d <- dat %>%
    filter(
      (grpA == gA & tmA == tA & grpB == gB & tmB == tB) |
        (grpA == gB & tmA == tB & grpB == gA & tmB == tA)
    ) %>%
    {
      if (length(ex)) filter(., !(subjA %in% ex | subjB %in% ex)) else .
    }

  # 2) count matches in db; return empty template if none
  n_matches <- d %>%
    summarise(n = dplyr::n()) %>%
    collect() %>%
    pull(n)
  if (is.na(n_matches) || n_matches == 0L) {
    return(tmpl)
  }

  # 3) pair by subject or by dyad
  if (by == "subject") {
    d <- d %>%
      filter(subjA == subjB) %>%
      select(similarity, id = subjA)
  } else { # dyad
    d <- d %>%
      filter(dyadA == dyadB) %>%
      select(similarity, id = dyadA)
  }

  # 4) collect and return labelled tibble
  d <- d %>% collect()
  if (!nrow(d)) {
    return(tmpl)
  }

  tibble(label = label, block = block, id = d$id, similarity = d$similarity)
}

# build long df for the fixed cmp_specs
stab_long <- pmap_dfr(cmp_specs, pull_cmp_safe)

# ---------- 4) adults (predicts / controls) ----------
# CAVE: 64GB of RAM required at least
withr::with_preserve_seed({
  ps1 <- phip_convert(
    data_long_path = "data/PREDICTS_counts.parquet",
    backend = "duckdb",
    peptide_library = TRUE,
    subject_id = "SUBJECTID",
    sample_id = "sample_id",
    peptide_id = "peptide_name",
    fold_change = "fold_change",
    timepoint = "timepoint",
    materialise_table = TRUE,
    auto_expand = TRUE,
    n_cores = 12
  )
})

ps1 %<>% filter(GROUP == "Healthy controls")

# those are old phiper utils i used earlier
source("scripts/99-utils.R")

stab_df2 <- compute_repertoire_stability(
  x = ps1, group_col = NULL, time_col = "timepoint",
  similarity = "kulczynski", drop_self = TRUE
)

# map adults time --> labels (keep as before)
time_map <- c(`-4` = "6 years", `-2` = "8 years", `0` = "10 years")

adults_long <- stab_df2 %>%
  mutate(label = recode(as.character(time), !!!time_map)) %>%
  filter(!is.na(label)) %>%
  transmute(label, block = "Adults", id = subject_id, similarity)

# ---------- 5) join and axis ordering ----------
base_order <- label_order
label_order2 <- c(base_order, unname(time_map[c("-4", "-2", "0")]))

stab_all <- bind_rows(stab_long, adults_long) %>%
  mutate(
    label = factor(label, levels = label_order2),
    block = factor(block, levels = c(
      "Pregnancy", "Infants",
      "Breast milk", "Adults"
    ))
  )

# ---------- 6) colors by label (left->right) ----------
raw_cols10 <- c(
  "#f54933", "#f5335d", "#f57f98", "#ffff00",
  "#05f5a5", "#33f5d2", "#1c6c4c",
  "#cae1ff", "#00bfff", "#1a1a8a"
)

# fill missing three middle colours (breast milk-ish)
auto_fill <- c("#e0087f", "#d037de", "#9a37de")

cols_full <- c(raw_cols10, auto_fill)
stopifnot(length(cols_full) == length(label_order2))

label_pal <- setNames(cols_full, label_order2)

# ---------- 7) main figure (phip style) ----------
adult_start <- which(levels(stab_all$label) %in%
  unname(time_map[c("-4", "-2", "0")]))[1]
vline_x <- adult_start - 0.5

p <- ggplot(stab_all, aes(
  x = label, y = similarity,
  fill = label, color = label
)) +
  geom_boxplot(
    width = 0.70, outlier.shape = NA, linewidth = 0.3,
    alpha = 0.85, color = "black"
  ) +
  geom_point(
    position = position_jitter(width = 0.24, height = 0, seed = 1),
    size = 0.05, alpha = 0.3, color = "black", inherit.aes = FALSE,
    data = stab_all, mapping = aes(x = label, y = similarity)
  ) +
  scale_fill_manual(values = label_pal, guide = "none") +
  scale_color_manual(values = label_pal, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Kulczynski similarity") +
  theme_phip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 20, 10, 10),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20)
  ) +
  geom_vline(xintercept = vline_x, linetype = "dotted") +
  scale_y_continuous(breaks = seq(0, 1, 0.1))
p

# helper to add header bands
add_header <- function(p, lab, x0, x1, col) {
  p +
    annotate("segment",
      x = x0, xend = x1, y = 1.04, yend = 1.04,
      linewidth = 3, colour = col
    ) +
    annotate("label",
      x = (x0 + x1) / 2, y = 1.07, label = lab,
      vjust = 0, size = 6, fontface = "bold",
      fill = "white", label.size = 0
    )
}

rng <- list(
  Pregnancy     = c(0.5, 4),
  Infants       = c(4, 7.5),
  `Breast milk` = c(7.5, 10.5),
  Adults        = c(vline_x, length(levels(stab_all$label)) + 0.5)
)

# expand plotting area upward to show headers
p <- p + coord_cartesian(ylim = c(0, 1.10), clip = "off") +
  theme(plot.margin = margin(10, 10, 10, 10))

p <- p |>
  add_header("Pregnancy", rng$Pregnancy[1], rng$Pregnancy[2], "#f54933") |>
  add_header("Infants", rng$Infants[1], rng$Infants[2], "#05f5a5") |>
  add_header(
    "Breast milk", rng$`Breast milk`[1],
    rng$`Breast milk`[2], "#0072B2"
  ) |>
  add_header("Adults", rng$Adults[1], rng$Adults[2], "#c46be8")

ggsave("results/stab-rich_plots/stability_boxplot.png", p,
  height = 8, width = 11, units = "cm", dpi = 300, bg = "white"
)

# ---------- 3b) collect available timepoints per group (db safe) --------------
times_by_group <- norm_pairs %>%
  select(group = grpA, time = tmA) %>%
  union_all(norm_pairs %>% select(group = grpB, time = tmB)) %>%
  distinct() %>%
  collect() %>%
  arrange(group, time) %>%
  group_by(group) %>%
  summarise(times = list(sort(unique(time))), .groups = "drop")

# helper returning vector of times for a group
get_times <- function(g) {
  i <- which(times_by_group$group == g)
  if (length(i)) times_by_group$times[[i[1]]] else character()
}

bm_times <- get_times("mom_milk")
kid_times <- get_times("kid_serum")
mom_times <- get_times("mom_serum")

# ---------- 3c) ordering function: use mapped printed labels (m0 -> b) --------
time_levels_all <- c("P12", "P28", "B", "M1", "M3", "M6", "M12")
ord_time <- function(x) {
  factor(map_time_label(x),
    levels = intersect(
      time_levels_all,
      time_levels_all %in% unique(map_time_label(x))
    )
  )
}

# ---------- bm vs kid_serum (full grid) ----------
grid_BM_I <- tidyr::expand_grid(time_bm = bm_times, time_kid = kid_times)

cmp_specs_BM_I_tmp <- grid_BM_I %>%
  arrange(ord_time(time_bm), ord_time(time_kid)) %>%
  mutate(
    # use mapped printed labels for readability (m0 -> b)
    label = paste0(
      "BM_", map_time_label(time_bm), " vs I_",
      map_time_label(time_kid)
    ),
    group_a = "mom_milk", time_a = time_bm,
    group_b = "kid_serum", time_b = time_kid,
    match_by = "dyad", block = "BM vs kid_serum"
  )

cmp_specs_BM_I <- cmp_specs_BM_I_tmp %>%
  transmute(label, group_a, time_a, group_b, time_b, match_by, block)

ord_BM_I <- cmp_specs_BM_I_tmp$label

# ---------- bm vs mom_serum (full grid) ----------
grid_BM_M <- tidyr::expand_grid(time_bm = bm_times, time_mom = mom_times)

cmp_specs_BM_M_tmp <- grid_BM_M %>%
  arrange(ord_time(time_bm), ord_time(time_mom)) %>%
  mutate(
    label = paste0(
      "BM_", map_time_label(time_bm), " vs M_",
      map_time_label(time_mom)
    ),
    group_a = "mom_milk", time_a = time_bm,
    group_b = "mom_serum", time_b = time_mom,
    match_by = "subject", block = "BM vs mom_serum"
  )

cmp_specs_BM_M <- cmp_specs_BM_M_tmp %>%
  transmute(label, group_a, time_a, group_b, time_b, match_by, block)

ord_BM_M <- cmp_specs_BM_M_tmp$label

# -----------------------
# add: mom_serum vs kid_serum (full grid)
# -----------------------
grid_mom_kid <- tidyr::expand_grid(time_mom = mom_times, time_kid = kid_times)

cmp_specs_mom_kid_tmp <- grid_mom_kid %>%
  arrange(ord_time(time_mom), ord_time(time_kid)) %>%
  mutate(
    label = paste0(
      "M_", map_time_label(time_mom), " vs I_",
      map_time_label(time_kid)
    ),
    group_a = "mom_serum", time_a = time_mom,
    group_b = "kid_serum", time_b = time_kid,
    match_by = "dyad", # match mom <-> kid by dyad
    block = "mom_serum vs kid_serum"
  )

cmp_specs_mom_kid <- cmp_specs_mom_kid_tmp %>%
  transmute(label, group_a, time_a, group_b, time_b, match_by, block)

ord_mom_kid <- cmp_specs_mom_kid_tmp$label

# ---------- fetch pairs and assign factor ordering ----------
stab_BM_I <- purrr::pmap_dfr(cmp_specs_BM_I, pull_cmp_safe,
  dat = norm_pairs
) %>%
  mutate(label = factor(label, levels = ord_BM_I), block = factor(block))

stab_BM_M <- purrr::pmap_dfr(cmp_specs_BM_M, pull_cmp_safe,
  dat = norm_pairs
) %>%
  mutate(label = factor(label, levels = ord_BM_M), block = factor(block))

stab_mom_kid <- purrr::pmap_dfr(cmp_specs_mom_kid, pull_cmp_safe,
  dat = norm_pairs
) %>%
  mutate(label = factor(label, levels = ord_mom_kid), block = factor(block))

# ---------- 3d) three plots: bm vs kid, bm vs mom, mom vs kid ----------
col_kid <- "#05f5a5"
col_mom <- "#f54933"
col_cross <- "orange"

p_BM_I <- ggplot(stab_BM_I, aes(x = label, y = similarity)) +
  geom_boxplot(
    fill = col_kid, color = "black", width = 0.70,
    outlier.shape = NA, linewidth = 0.3, alpha = 0.9
  ) +
  geom_point(
    data = stab_BM_I, aes(x = label, y = similarity),
    position = position_jitter(width = 0.24, height = 0, seed = 1),
    size = 0.05, alpha = 0.30, color = "black", inherit.aes = FALSE
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  labs(
    x = NULL, y = "Kulczynski similarity",
    title = "Breast milk vs kid serum (all timepoint pairs;
       matching by dyad)"
  ) +
  theme_phip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 20, 10, 10),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20)
  )

p_BM_M <- ggplot(stab_BM_M, aes(x = label, y = similarity)) +
  geom_boxplot(
    fill = col_mom, color = "black", width = 0.70,
    outlier.shape = NA, linewidth = 0.3, alpha = 0.9
  ) +
  geom_point(
    data = stab_BM_M, aes(x = label, y = similarity),
    position = position_jitter(width = 0.24, height = 0, seed = 1),
    size = 0.05, alpha = 0.30, color = "black", inherit.aes = FALSE
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  labs(
    x = NULL, y = "Kulczynski similarity",
    title = "Breast milk vs mom serum (all timepoint pairs;
       matching by subject)"
  ) +
  theme_phip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 20, 10, 10),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20)
  )

p_mom_kid <- ggplot(stab_mom_kid, aes(x = label, y = similarity)) +
  geom_boxplot(
    fill = "#ffb84d", color = "black", width = 0.70,
    outlier.shape = NA, linewidth = 0.3, alpha = 0.9
  ) +
  geom_point(
    data = stab_mom_kid, aes(x = label, y = similarity),
    position = position_jitter(width = 0.24, height = 0, seed = 1),
    size = 0.05, alpha = 0.30, color = "black", inherit.aes = FALSE
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  labs(
    x = NULL, y = "Kulczynski similarity",
    title = "Mom serum vs kid serum (all timepoint pairs;
       matching by dyad)"
  ) +
  theme_phip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 20, 10, 10),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16)
  )

# preview plots
p_BM_I
p_BM_M
p_mom_kid

# save the three plots
ggsave("results/stab-rich_plots/stability_BM_vs_kid_serum.png", p_BM_I,
  height = 8, width = 10, units = "cm", dpi = 300, bg = "white"
)
ggsave("results/stab-rich_plots/stability_BM_vs_mom_serum.png", p_BM_M,
  height = 8, width = 10, units = "cm", dpi = 300, bg = "white"
)
ggsave("results/stab-rich_plots/stability_mom_vs_kid_serum.png", p_mom_kid,
  height = 8, width = 10, units = "cm", dpi = 300, bg = "white"
)

# --- BM M6 vs everything in mom_serum and kid_serum (single combined plot) ----
# 1) build grids: always bm_m6 vs each time in mom and kid groups
grid_BM6_M <- tibble::tibble(time_bm = "M6") %>%
  tidyr::crossing(time_mom = mom_times)
grid_BM6_I <- tibble::tibble(time_bm = "M6") %>%
  tidyr::crossing(time_kid = kid_times)

# 2) temporary specs with mapped label strings for bm6
cmp_BM6_M_tmp <- grid_BM6_M %>%
  arrange(ord_time(time_bm), ord_time(time_mom)) %>%
  mutate(
    label = paste0("BM_M6 vs M_", map_time_label(time_mom)),
    group_a = "mom_milk", time_a = time_bm,
    group_b = "mom_serum", time_b = time_mom,
    match_by = "subject",
    block = "BM_M6 vs mom_serum",
    type = "mom_serum"
  )

cmp_BM6_I_tmp <- grid_BM6_I %>%
  arrange(ord_time(time_bm), ord_time(time_kid)) %>%
  mutate(
    label = paste0("BM_M6 vs I_", map_time_label(time_kid)),
    group_a = "mom_milk", time_a = time_bm,
    group_b = "kid_serum", time_b = time_kid,
    match_by = "dyad",
    block = "BM_M6 vs kid_serum",
    type = "kid_serum"
  )

# 3) build final extra block: mom_serum P28 vs all kid_serum times (dyad)
grid_momP28_kid <- tibble::tibble(time_mom = "P28") %>%
  tidyr::crossing(time_kid = kid_times)

cmp_momP28_kid_tmp <- grid_momP28_kid %>%
  arrange(ord_time(time_mom), ord_time(time_kid)) %>%
  mutate(
    label = paste0(
      "M_", map_time_label(time_mom), " vs I_",
      map_time_label(time_kid)
    ),
    group_a = "mom_serum", time_a = time_mom,
    group_b = "kid_serum", time_b = time_kid,
    match_by = "dyad",
    block = "M_P28 vs kid_serum",
    type = "momP28_vs_kid"
  )

# 5) final specs for pull_cmp_safe (drop 'type' column)
cmp_BM6_M <- cmp_BM6_M_tmp %>%
  transmute(label, group_a, time_a, group_b, time_b, match_by, block)

cmp_BM6_I <- cmp_BM6_I_tmp %>%
  transmute(label, group_a, time_a, group_b, time_b, match_by, block)

cmp_momP28_kid <- cmp_momP28_kid_tmp %>%
  transmute(label, group_a, time_a, group_b, time_b, match_by, block)

# 6) fetch pairs from db
stab_BM6_M <- purrr::pmap_dfr(cmp_BM6_M, pull_cmp_safe, dat = norm_pairs)
stab_BM6_I <- purrr::pmap_dfr(cmp_BM6_I, pull_cmp_safe, dat = norm_pairs)
stab_momP28_kid <- purrr::pmap_dfr(cmp_momP28_kid, pull_cmp_safe,
  dat = norm_pairs
)

# 7) stitch together and set explicit ordering + colors
# note: we pick the exact order you requested:
#   mom block:  M_P12, M_P28, M_B
#   kid block:  I_B, I_M3, I_M12
#   momP28->kid block: same kid order (3 orange shades)

# stitch data
stab_bm6_all <- bind_rows(stab_BM6_M, stab_BM6_I, stab_momP28_kid) %>%
  left_join(bind_rows(
    cmp_BM6_M_tmp %>% select(label, type),
    cmp_BM6_I_tmp %>% select(label, type),
    cmp_momP28_kid_tmp %>% select(label, type)
  ), by = "label")

# define desired time orders (map_time_label must be available)
mom_time_order <- c("P12", "P28", "B")
kid_time_order <- c("B", "M3", "M12")

# extract labels in desired order (defensive: intersect)
mom_labels <- cmp_BM6_M_tmp %>%
  filter(map_time_label(time_mom) %in% mom_time_order) %>%
  mutate(ord = match(map_time_label(time_mom), mom_time_order)) %>%
  arrange(ord) %>%
  pull(label)

kid_labels <- cmp_BM6_I_tmp %>%
  filter(map_time_label(time_kid) %in% kid_time_order) %>%
  mutate(ord = match(map_time_label(time_kid), kid_time_order)) %>%
  arrange(ord) %>%
  pull(label)

momP28_labels <- cmp_momP28_kid_tmp %>%
  filter(map_time_label(time_kid) %in% kid_time_order) %>%
  mutate(ord = match(map_time_label(time_kid), kid_time_order)) %>%
  arrange(ord) %>%
  pull(label)

# final x-axis order: mom block, kid block, then momP28->kid block
ord_labels <- c(mom_labels, kid_labels, momP28_labels)

# force factor levels to that order
stab_bm6_all <- stab_bm6_all %>%
  mutate(
    label = factor(label, levels = ord_labels),
    type  = factor(type, levels = c("mom_serum", "kid_serum", "momP28_vs_kid"))
  )

# 8) build color map consistent with requested order
cols_mom_by_time <- c(P12 = "#f54933", P28 = "#f5335d", B = "#f57f98")
cols_kid_by_time <- c(B = "#05f5a5", M3 = "#33f5d2", M12 = "#1c6c4c")

# three orange shades for the momP28->kid block (light -> dark)
oranges <- c("#ffcf99", "#ff9f33", "#ff6600")

# map colors to labels (keep explicit order)
lab_col_mom <- setNames(
  cols_mom_by_time[map_time_label(c("P12", "P28", "B"))],
  paste0(
    "BM_M6 vs M_",
    map_time_label(c("P12", "P28", "B"))
  )
)
lab_col_kid <- setNames(
  cols_kid_by_time[map_time_label(c("B", "M3", "M12"))],
  paste0(
    "BM_M6 vs I_",
    map_time_label(c("B", "M3", "M12"))
  )
)
lab_col_momP28 <- setNames(
  oranges,
  paste0(
    "M_", map_time_label("P28"), " vs I_",
    map_time_label(c("B", "M3", "M12"))
  )
)

# combine and ensure names match ord_labels; unknown labels get a safe fallback
label_colors <- c(lab_col_mom, lab_col_kid, lab_col_momP28)
missing_lbls <- setdiff(ord_labels, names(label_colors))
if (length(missing_lbls)) {
  label_colors[missing_lbls] <- "#bbbbbb"
}
# reorder label_colors to match ord_labels
label_colors <- label_colors[ord_labels]

# 9) build plot (uses the explicit ordering and colors above)
p_bm6 <- ggplot(stab_bm6_all, aes(x = label, y = similarity, fill = label)) +
  geom_boxplot(
    width = 0.70, outlier.shape = NA, linewidth = 0.3, alpha = 0.9,
    color = "black"
  ) +
  geom_point(
    position = position_jitter(width = 0.24, height = 0, seed = 1),
    size = 0.05, alpha = 0.30, color = "black"
  ) +
  scale_fill_manual(values = label_colors, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  labs(x = NULL, y = "Kulczynski similarity") +
  theme_phip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 20, 10, 10),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20)
  )

# 10) add band headers (three sections) and save
n_mom <- length(mom_labels)
n_kid <- length(kid_labels)
n_momP28 <- length(momP28_labels)

p_bm6 <- p_bm6 +
  coord_cartesian(ylim = c(0, 1.10), clip = "off") +
  theme(plot.margin = margin(10, 10, 10, 10))

p_bm6 <- p_bm6 |>
  add_header("BM_M6 vs mom_serum", 0.5, n_mom + 0.5, col_mom) |>
  add_header(
    "BM_M6 vs kid_serum", n_mom + 0.5, n_mom + n_kid + 0.5,
    col_kid
  ) |>
  add_header(
    "M_P28 vs kid_serum", n_mom + n_kid + 0.5,
    n_mom + n_kid + n_momP28 + 0.5, "#ff8000"
  )

# preview + save
p_bm6
ggsave("results/stab-rich_plots/stability_BM_M6_vs_all_ordered.png", p_bm6,
  height = 8, width = 10, units = "cm", dpi = 300, bg = "white"
)


# --- drop other temp tables (keep only the new data_long) ---
con <- ps_merged_box_bin$meta$con
all_tables <- DBI::dbListTables(con)

orig_tbl <- if (is.character(ps_merged_box_bin$data_long)) {
  ps_merged_box_bin$data_long
} else {
  dbplyr::remote_name(ps_merged_box_bin$data_long)
}

to_keep <- orig_tbl
to_drop <- setdiff(all_tables, to_keep)
if (length(to_drop) > 0) {
  for (t in to_drop) {
    DBI::dbExecute(con, paste0(
      "DROP TABLE IF EXISTS ",
      DBI::dbQuoteIdentifier(con, t)
    ))
  }
}

# refresh R reference, ramove objects and free memory
ps_merged_box_bin$data_long <- dplyr::tbl(con, orig_tbl)
rm(
  list = setdiff(
    ls(envir = .GlobalEnv, all.names = TRUE),
    c("ps_merged_box_bin", "seed_before")
  ),
  envir = .GlobalEnv
)

gc()
gc()
