long_to_wide_exist <- function(
    x,
    row_cols   = c("subject_id","big_group","timepoint_factor"),
    col_col    = "peptide_id",
    value_col  = "exist",
    sep        = "__",
    sparse     = TRUE
) {
  # if it's a dbplyr/duckdb lazy table, pull only the needed columns first
  if (inherits(x, "tbl_sql")) {
    x <- dplyr::select(x, dplyr::all_of(c(row_cols, col_col, value_col))) %>%
      dplyr::collect()
  }

  if (rlang::is_installed("data.table")) {
    DT <- data.table::as.data.table(x)

    # build row_id once; keep only needed cols
    data.table::set(DT, j = "row_id",
                    value = do.call(paste, c(DT[, ..row_cols], list(sep = sep))))
    DT <- DT[, c("row_id", col_col, value_col), with = FALSE]

    # force integer 0/1
    data.table::set(DT, j = value_col,
                    value = as.integer(DT[[value_col]] > 0))

    if (isTRUE(sparse)) {
      # sparse path (recommended for many peptides)
      # duplicates are *summed* by sparseMatrix; threshold back to 0/1
      row_f <- factor(DT$row_id)
      col_f <- factor(DT[[col_col]])
      M <- Matrix::sparseMatrix(
        i = as.integer(row_f),
        j = as.integer(col_f),
        x = DT[[value_col]],
        dims = c(nlevels(row_f), nlevels(col_f)),
        dimnames = list(levels(row_f), levels(col_f))
      )
      return((M > 0) + 0L) # dgCMatrix with 0/1
    } else {
      # dense matrix via dcast (watch memory!)
      form <- reformulate(col_col, response = "row_id")
      wide_dt <- data.table::dcast(
        DT, formula = form,
        value.var = value_col,
        fun.aggregate = max,  # collapse duplicates to 0/1
        fill = 0L
      )
      rn <- wide_dt[["row_id"]]
      mat <- as.matrix(wide_dt[, setdiff(names(wide_dt), "row_id"), with = FALSE])
      rownames(mat) <- rn
      storage.mode(mat) <- "integer"
      return(mat)
    }
  } else {
    # tidyr fallback (will be slower and memory-heavier)
    wide <- x |>
      dplyr::mutate(
        row_id = do.call(paste, c(dplyr::pick(dplyr::all_of(row_cols)), list(sep = sep))),
        !!value_col := as.integer(.data[[value_col]] > 0)
      ) |>
      dplyr::select(row_id, dplyr::all_of(c(col_col, value_col))) |>
      dplyr::group_by(row_id, .data[[col_col]]) |>
      dplyr::summarise(!!value_col := max(.data[[value_col]]), .groups = "drop") |>
      tidyr::pivot_wider(
        id_cols     = row_id,
        names_from  = dplyr::all_of(col_col),
        values_from = dplyr::all_of(value_col),
        values_fill = 0L
      )

    rn <- wide$row_id
    mat <- as.matrix(wide[, setdiff(names(wide), "row_id"), drop = FALSE])
    rownames(mat) <- rn
    storage.mode(mat) <- "integer"
    if (isTRUE(sparse)) {
      return(Matrix::Matrix(mat, sparse = TRUE))
    } else {
      return(mat)
    }
  }
}

kulczynski_dist <- function(M) {
  # ensure sparse binary
  M <- drop0(M)
  M <- as(M != 0, "dgCMatrix")

  # intersections (dense n x n; fine for n = nrow(M))
  A   <- as.matrix(tcrossprod(M))     # a_ij
  r   <- Matrix::rowSums(M)           # r_i
  r_s <- ifelse(r == 0, 1, r)         # safe denom to avoid /0

  # 0.5 * (a/r_i + a/r_j)
  K <- 0.5 * (sweep(A, 1, r_s, "/") + sweep(A, 2, r_s, "/"))

  # distances
  D <- 1 - K
  diag(D) <- 0
  D
}

# ==============================================================================
# Repertoire stability: compute (Jaccard / Sorensen / Kulczynski / Bray-Curtis / Phi) --> tibble
# ==============================================================================

#' @title Compute repertoire stability over time (Jaccard, Sorensen, Kulczynski, Bray-Curtis, Phi)
#'
#' @description For each **subject**, within each **group**, compute the
#' similarity between the set of *present* peptides at each time point and that
#' subject’s **baseline** (their earliest time point). Supported:
#' - **"jaccard"**    : \eqn{|A ∩ B| / |A ∪ B|}
#' - **"sorensen"**   : \eqn{2|A ∩ B| / (|A| + |B|)}  (Dice)
#' - **"kulczynski"** : \eqn{0.5 * [ |A∩B|/(|A∩B|+|A\B|) + |A∩B|/(|A∩B|+|B\A|) ]}
#' - **"bray_curtis"** (binary here): identical to Sørensen/Dice
#' - **"phi"** / **"pearson"** : phi (Pearson corr.) for two binarne wektory; zakres [-1, 1]
#'
#' Baseline self-comparisons (time == baseline) equal 1 by definition for set-based
#' similarities; if `drop_self = TRUE` these rows are removed, a marker is kept.
#'
#' @param x A `<phip_data>` object. **Must be longitudinal**.
#' @param group_col Character; grouping variable in `x$data_long`.
#' @param time_col Character; numeric/Date/POSIX time variable.
#' @param similarity One of: "jaccard","sorensen","kulczynski","bray_curtis","phi","pearson".
#' @param drop_self Logical; drop baseline rows (default `TRUE`).
#' @param auto_expand Logical; if `TRUE`, require n0>=1 & nt>=1.
#'
#' @return tibble: `subject_id`, `group`, `time`, `similarity` (+ attributes).
#' @export
compute_repertoire_stability <- function(
    x,
    group_col,
    time_col,
    similarity  = c("jaccard", "sorensen", "kulczynski", "bray_curtis", "phi", "pearson"),
    drop_self   = TRUE,
    auto_expand = FALSE,  # require n0>=1 & nt>=1 when TRUE
    full_cross  = FALSE   # NEW: all pairwise (group,time) combos per subject
) {
  stopifnot(inherits(x, "phip_data"))
  similarity <- match.arg(similarity)

  group_col_local <- if (is.null(group_col)) ".phip__group" else group_col

  .ph_with_timing(
    headline = "Computing repertoire similarity (<phip_data>)",
    step     = sprintf("similarity: %s; mode: %s",
                       similarity, if (full_cross) "all-pairs" else "baseline"),
    expr = {
      if (!isTRUE(x$meta$longitudinal)) {
        .ph_abort(
          headline = "Longitudinal data required.",
          step     = "meta$longitudinal check",
          bullets  = "x$meta$longitudinal is FALSE"
        )
      }

      # --- wymogi kolumn ------------------------------------------------------
      need <- c("subject_id", "sample_id", "peptide_id", "exist", time_col)
      if (!is.null(group_col)) need <- c(need, group_col)
      miss <- setdiff(need, colnames(x$data_long))
      .chk_cond(length(miss) > 0,
                sprintf("Missing required columns in data_long: %s",
                        paste(miss, collapse = ", ")))

      .data <- rlang::.data
      tcol  <- rlang::sym(time_col)
      grp   <- rlang::sym(group_col_local)

      # --- presence-only widok (exist==1) + syntetyczna grupa "all" ----------
      present <- if (is.null(group_col)) {
        x$data_long |>
          dplyr::filter(.data$exist == 1L) |>
          dplyr::select(.data$subject_id, .data$sample_id, .data$peptide_id, !!tcol) |>
          dplyr::filter(!is.na(!!tcol)) |>
          dplyr::mutate(!!grp := "all")
      } else {
        x$data_long |>
          dplyr::filter(.data$exist == 1L) |>
          dplyr::select(.data$subject_id, .data$sample_id, .data$peptide_id, !!grp, !!tcol) |>
          dplyr::filter(!is.na(!!tcol))
      }

      # =======================================================================
      # BRANCH A: full_cross = TRUE  --> wszystkie pary (group,time)
      # =======================================================================
      if (isTRUE(full_cross)) {

        # zbiory: unikatowe peptydy per (subject, group, time)
        sets <- present |>
          dplyr::distinct(.data$subject_id, !!grp, !!tcol, .data$peptide_id, .data$sample_id)

        # rozmiary zbiorów
        sizes <- sets |>
          dplyr::count(.data$subject_id, !!grp, !!tcol, name = "n")

        # self-join po peptydzie w obrębie osoby; filtrujemy do par unikalnych (key1 < key2)
        sets_l <- sets |>
          dplyr::rename(group1 = !!grp, time1 = !!tcol,
                        sample_id1 = .data$sample_id, peptide_id1 = .data$peptide_id)
        sets_r <- sets |>
          dplyr::rename(group2 = !!grp, time2 = !!tcol,
                        sample_id2 = .data$sample_id, peptide_id2 = .data$peptide_id)

        inters <- sets_l |>
          dplyr::inner_join(
            sets_r,
            by = c("subject_id", "peptide_id1" = "peptide_id2")
          ) |>
          dplyr::filter(.data$group1 != .data$group2 | .data$time1 != .data$time2) |>
          dplyr::mutate(
            key1 = paste0(.data$group1, "||", .data$time1),
            key2 = paste0(.data$group2, "||", .data$time2)
          ) |>
          dplyr::filter(.data$key1 < .data$key2) |>
          dplyr::count(.data$subject_id, .data$group1, .data$time1,
                       .data$group2, .data$time2, name = "a")

        # dołącz rozmiary
        pairs <- inters |>
          dplyr::left_join(
            sizes |> dplyr::rename(group1 = !!grp, time1 = !!tcol, n1 = .data$n),
            by = c("subject_id","group1","time1")
          ) |>
          dplyr::left_join(
            sizes |> dplyr::rename(group2 = !!grp, time2 = !!tcol, n2 = .data$n),
            by = c("subject_id","group2","time2")
          ) |>
          dplyr::mutate(
            a = dplyr::coalesce(.data$a, 0L),
            n1 = dplyr::coalesce(.data$n1, 0L),
            n2 = dplyr::coalesce(.data$n2, 0L),
            b = pmax(.data$n1 - .data$a, 0L),
            c = pmax(.data$n2 - .data$a, 0L)
          )

        # phi potrzebuje N – policz rozmiary przestrzeni cech per próbka
        if (identical(similarity, "phi") || identical(similarity, "pearson")) {
          timesamples <- x$data_long |>
            dplyr::filter(!is.na(!!tcol)) |>
            dplyr::select(.data$subject_id, !!grp, !!tcol, .data$sample_id) |>
            dplyr::distinct()

          sample_sizes <- x$data_long |>
            dplyr::group_by(.data$sample_id) |>
            dplyr::summarise(N = dplyr::n_distinct(.data$peptide_id), .groups = "drop")

          pairs <- pairs |>
            # próbki odpowiadające parze (group1,time1) / (group2,time2)
            dplyr::left_join(
              timesamples |>
                dplyr::rename(group1 = !!grp, time1 = !!tcol, sample_id1 = .data$sample_id),
              by = c("subject_id","group1","time1")
            ) |>
            dplyr::left_join(
              timesamples |>
                dplyr::rename(group2 = !!grp, time2 = !!tcol, sample_id2 = .data$sample_id),
              by = c("subject_id","group2","time2")
            ) |>
            dplyr::left_join(sample_sizes, by = c("sample_id1" = "sample_id")) |>
            dplyr::rename(N1 = .data$N) |>
            dplyr::left_join(sample_sizes, by = c("sample_id2" = "sample_id")) |>
            dplyr::rename(N2 = .data$N) |>
            dplyr::mutate(
              N   = pmin(.data$N1, .data$N2, na.rm = TRUE),
              n00 = pmax(.data$N - (.data$a + .data$b + .data$c), 0),
              n1d = .data$a + .data$b,
              nd1 = .data$a + .data$c,
              n0d = .data$N - .data$n1d,
              nd0 = .data$N - .data$nd1,
              den = sqrt(pmax(.data$n1d,0) * pmax(.data$n0d,0) *
                           pmax(.data$nd1,0) * pmax(.data$nd0,0)),
              similarity = dplyr::if_else(.data$den > 0,
                                          (.data$a * .data$n00 - .data$b * .data$c) / .data$den,
                                          NA_real_)
            ) |>
            dplyr::select(-n00,-n1d,-nd1,-n0d,-nd0,-den,-N1,-N2)
          if (!isTRUE(x$meta$full_cross)) {
            .ph_warn("Phi/Pearson requested but meta$full_cross is FALSE; N may be underestimated.")
          }
        } else if (identical(similarity, "jaccard")) {
          pairs <- pairs |>
            dplyr::mutate(
              similarity = dplyr::if_else((.data$n1 + .data$n2 - .data$a) > 0,
                                          .data$a / (.data$n1 + .data$n2 - .data$a),
                                          NA_real_)
            )
        } else if (identical(similarity, "sorensen") || identical(similarity, "bray_curtis")) {
          pairs <- pairs |>
            dplyr::mutate(
              similarity = dplyr::if_else((.data$n1 + .data$n2) > 0,
                                          (2 * .data$a) / (.data$n1 + .data$n2),
                                          NA_real_)
            )
        } else if (identical(similarity, "kulczynski")) {
          pairs <- pairs |>
            dplyr::mutate(
              term1 = dplyr::if_else((.data$a + .data$b) > 0, .data$a / (.data$a + .data$b), NA_real_),
              term2 = dplyr::if_else((.data$a + .data$c) > 0, .data$a / (.data$a + .data$c), NA_real_),
              similarity = dplyr::coalesce(0.5 * (.data$term1 + .data$term2),
                                           dplyr::if_else((.data$a + .data$b + .data$c) == 0, 1.0, 0.0))
            ) |>
            dplyr::select(-term1, -term2)
        }

        out <- pairs |>
          dplyr::select(.data$subject_id, .data$group1, .data$time1, .data$group2, .data$time2, .data$similarity) |>
          dplyr::collect()
        class(out) <- c("phip_repertoire_pairs", class(out))
        return(out)
      }

      # =======================================================================
      # BRANCH B: baseline (oryginalny tryb)
      # =======================================================================

      # ---- baseline (per subject within group) & sets ------------------------
      baseline_tbl <- present |>
        dplyr::group_by(.data$subject_id, !!grp) |>
        dplyr::summarise(baseline_time = min(!!tcol, na.rm = TRUE), .groups = "drop")

      base_set <- present |>
        dplyr::inner_join(baseline_tbl, by = c("subject_id", group_col_local)) |>
        dplyr::filter(!!tcol == .data$baseline_time) |>
        dplyr::distinct(.data$subject_id, !!grp, .data$peptide_id)

      tp_set <- present |>
        dplyr::distinct(.data$subject_id, !!grp, !!tcol, .data$peptide_id)

      bsize <- base_set |>
        dplyr::count(.data$subject_id, !!grp, name = "n0")
      tsize <- tp_set |>
        dplyr::count(.data$subject_id, !!grp, !!tcol, name = "nt")
      ints  <- tp_set |>
        dplyr::inner_join(base_set,
                          by = c("subject_id", group_col_local, "peptide_id")) |>
        dplyr::count(.data$subject_id, !!grp, !!tcol, name = "n_int")

      stab <- tsize |>
        dplyr::left_join(ints,  by = c("subject_id", group_col_local, time_col)) |>
        dplyr::left_join(bsize, by = c("subject_id", group_col_local)) |>
        dplyr::mutate(
          n_int = dplyr::coalesce(.data$n_int, 0L),
          nt    = dplyr::coalesce(.data$nt,    0L),
          n0    = dplyr::coalesce(.data$n0,    0L)
        ) |>
        dplyr::left_join(baseline_tbl, by = c("subject_id", group_col_local))

      if (isTRUE(auto_expand)) {
        n_before <- nrow(stab)
        stab <- stab |> dplyr::filter(.data$n0 >= 1L, .data$nt >= 1L)
        .ph_log_info("Auto-expand filter (require n0 >= 1 & nt >= 1)",
                     bullets = sprintf("removed rows: %d", n_before - nrow(stab)))
      }

      # ---- similarity (baseline) --------------------------------------------
      if (identical(similarity, "jaccard")) {
        stab <- stab |>
          dplyr::mutate(similarity = pmin(
            pmax(.data$n_int / pmax(.data$n0 + .data$nt - .data$n_int, 1L), 0), 1
          ))

      } else if (identical(similarity, "sorensen")) {
        stab <- stab |>
          dplyr::mutate(similarity = pmin(
            pmax((2 * .data$n_int) / pmax(.data$n0 + .data$nt, 1L), 0), 1
          ))

      } else if (identical(similarity, "bray_curtis")) {
        stab <- stab |>
          dplyr::mutate(similarity = pmin(
            pmax((2 * .data$n_int) / pmax(.data$n0 + .data$nt, 1L), 0), 1
          ))

      } else if (identical(similarity, "kulczynski")) {
        stab <- stab |>
          dplyr::mutate(
            a = .data$n_int,
            b = pmax(.data$n0 - .data$n_int, 0L),
            c = pmax(.data$nt - .data$n_int, 0L),
            term1 = dplyr::if_else((a + b) > 0, a / (a + b), NA_real_),
            term2 = dplyr::if_else((a + c) > 0, a / (a + c), NA_real_),
            similarity = dplyr::coalesce(0.5 * (term1 + term2),
                                         dplyr::if_else((a + b + c) == 0, 1.0, 0.0))
          ) |>
          dplyr::select(-a, -b, -c, -term1, -term2)

      } else { # phi / pearson
        timesamples <- x$data_long |>
          dplyr::filter(!is.na(!!tcol)) |>
          dplyr::select(.data$subject_id, !!grp, !!tcol, .data$sample_id) |>
          dplyr::distinct()

        baseline_sids <- timesamples |>
          dplyr::inner_join(
            baseline_tbl,
            by = c("subject_id", group_col_local, time_col = "baseline_time")
          ) |>
          dplyr::rename(baseline_sample_id = .data$sample_id) |>
          dplyr::select(.data$subject_id, !!grp, .data$baseline_sample_id)

        pair_sids <- timesamples |>
          dplyr::inner_join(baseline_sids, by = c("subject_id", group_col_local)) |>
          dplyr::rename(tp_sample_id = .data$sample_id)

        sample_sizes <- x$data_long |>
          dplyr::group_by(.data$sample_id) |>
          dplyr::summarise(N = dplyr::n_distinct(.data$peptide_id), .groups = "drop")

        N_map <- pair_sids |>
          dplyr::left_join(sample_sizes, by = c("tp_sample_id" = "sample_id")) |>
          dplyr::rename(N_tp = .data$N) |>
          dplyr::left_join(sample_sizes, by = c("baseline_sample_id" = "sample_id")) |>
          dplyr::rename(N_base = .data$N) |>
          dplyr::mutate(N = pmin(.data$N_tp, .data$N_base, na.rm = TRUE)) |>
          dplyr::select(.data$subject_id, !!grp, !!tcol, .data$N)

        stab <- stab |>
          dplyr::left_join(N_map, by = c("subject_id", group_col_local, time_col)) |>
          dplyr::mutate(
            a    = .data$n_int,
            b    = pmax(.data$n0 - .data$n_int, 0L),
            c    = pmax(.data$nt - .data$n_int, 0L),
            n00  = pmax(.data$N - (.data$a + .data$b + .data$c), 0),
            n1d  = .data$a + .data$b,
            nd1  = .data$a + .data$c,
            n0d  = .data$N - .data$n1d,
            nd0  = .data$N - .data$nd1,
            den  = sqrt(pmax(.data$n1d, 0) * pmax(.data$n0d, 0) *
                          pmax(.data$nd1, 0) * pmax(.data$nd0, 0)),
            similarity = dplyr::if_else(.data$den > 0,
                                        (.data$a * .data$n00 - .data$b * .data$c) / .data$den,
                                        NA_real_)
          ) |>
          dplyr::select(-a, -b, -c, -n00, -n1d, -nd1, -n0d, -nd0, -den)

        if (!isTRUE(x$meta$full_cross)) {
          .ph_warn("Phi/Pearson requested but meta$full_cross is FALSE; N may be underestimated. Similarity may be NA for some rows.")
        }
        if (isTRUE(any(is.na(stab$N)))) {
          n_na <- sum(is.na(stab$N))
          .ph_warn(sprintf("Phi/Pearson: missing N for %d row(s); similarity set to NA.", n_na))
        }
      }

      # --- finalize (baseline mode) ------------------------------------------
      out <- stab |>
        dplyr::select(.data$subject_id, !!grp, !!tcol, .data$similarity, .data$baseline_time) |>
        dplyr::collect() |>
        dplyr::rename(group = !!grp, time = !!tcol)

      base_per_group <- out |>
        dplyr::group_by(.data$group) |>
        dplyr::summarise(baseline_time = min(.data$baseline_time, na.rm = TRUE), .groups = "drop")

      if (isTRUE(drop_self)) {
        n_before <- nrow(out)
        out <- dplyr::filter(out, .data$time != .data$baseline_time)
        .ph_log_info("Dropping baseline self-comparisons",
                     bullets = sprintf("removed rows: %d", n_before - nrow(out)))
      }

      # baseline marker heights
      marker_tbl <- NULL
      if (nrow(out)) {
        marker_tbl <- out |>
          dplyr::select(-baseline_time) |>
          dplyr::inner_join(base_per_group, by = "group") |>
          dplyr::group_by(.data$group, .data$baseline_time) |>
          dplyr::summarise(
            next_time = {
              v <- .data$time[.data$time > .data$baseline_time]
              if (length(v)) min(v, na.rm = TRUE) else NA_real_
            },
            .groups = "drop"
          ) |>
          dplyr::left_join(
            out |> dplyr::group_by(.data$group, .data$time) |>
              dplyr::summarise(y_marker = mean(.data$similarity, na.rm = TRUE), .groups = "drop"),
            by = c("group", "next_time" = "time")
          ) |>
          dplyr::mutate(y_marker = dplyr::coalesce(.data$y_marker, 1)) |>
          dplyr::rename(x_marker = baseline_time)
      }

      class(out) <- c("phip_repertoire_stability", class(out))
      attr(out, "baseline_times")  <- base_per_group
      attr(out, "baseline_marker") <- marker_tbl
      attr(out, "drop_self")       <- isTRUE(drop_self)
      out
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ==============================================================================
# Repertoire stability: plot (takes tibble from compute_*)
# ==============================================================================

#' Plot repertoire stability over time (Jaccard or Sørensen)
#'
#' @description
#' Plot the output of [compute_repertoire_stability()] as group-wise curves
#' over time with confidence bands. If `drop_self = TRUE` was used at compute
#' time, the plot adds **one baseline marker per group**: a larger point at the
#' group baseline time (x), with y positioned at the **mean similarity of the
#' earliest non-baseline time** in that group, filled by group colour and with a
#' black border.
#'
#' @param stab_df A tibble returned by [compute_repertoire_stability()].
#'   Must contain columns `group`, `time`, `similarity`.
#' @param custom_colors Optional named vector of colors for groups.
#'   Defaults to PHIP palette.
#' @param continuous_mode `"gam"` (default), `"binned"`, or `"loess"`.
#' @param gam_k Requested GAM basis size per series (auto-shrinks to safe k).
#' @param ci_method `"model"` (default) or `"bootstrap"`.
#' @param ci_level Confidence level (default `0.95`).
#' @param boot_R Number of bootstrap replicates (default `500`).
#' @param boot_seed Optional integer seed for reproducibility.
#' @param boot_progress Logical; show textual progress bar during bootstrap
#'   (default `TRUE`).
#' @param ci_fill Ribbon fill color (default `"grey70"`).
#' @param ci_alpha Ribbon alpha (default `0.15`).
#' @param point_alpha Alpha for raw points overlay (default `0.25`; set `0` to hide).
#'
#' @return A `ggplot` object.
#' @export
plot_repertoire_stability <- function(
    stab_df,
    custom_colors        = NULL,
    # smoothing / ribbons
    continuous_mode      = c("gam","binned","loess"),
    gam_k                = 7,
    ci_method            = c("model","bootstrap"),
    ci_level             = 0.95,
    boot_R               = 500,
    boot_seed            = NULL,
    boot_progress        = TRUE,
    # visuals
    ci_fill              = "grey70",
    ci_alpha             = 0.15,
    point_alpha          = 0.25,
    # NEW: override fill colors for big baseline markers (named by group)
    marker_fill_colors   = NULL
) {
  stopifnot(is.data.frame(stab_df))
  continuous_mode <- match.arg(continuous_mode)
  ci_method       <- match.arg(ci_method)

  .ph_with_timing(
    headline = "Plotting repertoire stability (precomputed)",
    step     = sprintf("mode: %s; ci: %s", continuous_mode, ci_method),
    expr = {
      need <- c("group", "time", "similarity")
      miss <- setdiff(need, names(stab_df))
      .chk_cond(length(miss) > 0,
                sprintf("Missing required columns in stability tibble: %s",
                        paste(miss, collapse = ", ")))

      df <- tibble::as_tibble(stab_df)
      marker_tbl <- attr(stab_df, "baseline_marker", exact = TRUE)
      drop_self  <- isTRUE(attr(stab_df, "drop_self",     exact = TRUE))

      # lock group levels
      grp_levels <- if (!is.null(custom_colors) && length(names(custom_colors))) {
        names(custom_colors)
      } else unique(df$group)
      df$group <- factor(df$group, levels = grp_levels)
      if (is.data.frame(marker_tbl) && nrow(marker_tbl)) {
        marker_tbl$group <- factor(marker_tbl$group, levels = grp_levels)
      }

      # smooths per group
      split_list <- split(df, df$group)
      .ph_log_info("Fitting smooths per group",
                   bullets = c(sprintf("k requested: %d", gam_k),
                               sprintf("groups: %d", length(split_list))))
      preds <- purrr::map_dfr(split_list, function(d) {
        d <- as.data.frame(d)
        if (ci_method == "model") {
          out <- .gam_band_one(d, xvar = "time", yvar = "similarity",
                               k_req = gam_k, level = ci_level, nonneg = TRUE)
        } else {
          out <- .bootstrap_gam_one(d, xvar = "time", yvar = "similarity",
                                    k_req = gam_k, R = boot_R, level = ci_level,
                                    seed = boot_seed, progress = boot_progress,
                                    nonneg = TRUE)
        }
        out$group <- unique(d$group)
        out
      })
      preds$group <- factor(preds$group, levels = grp_levels)
      preds$.grp  <- preds$group

      # base plot (phiper style)
      p <- ggplot2::ggplot()
      if (point_alpha > 0) {
        p <- p +
          ggplot2::geom_point(
            data = df,
            ggplot2::aes(x = .data$time, y = .data$similarity, color = .data$group),
            alpha = point_alpha, size = 1.6, show.legend = TRUE
          )
      }
      p <- p +
        ggplot2::geom_ribbon(
          data = preds,
          ggplot2::aes(x = .data$.x, ymin = .data$lwr, ymax = .data$upr, group = .data$.grp),
          fill = ci_fill, alpha = ci_alpha, colour = NA, inherit.aes = FALSE
        ) +
        ggplot2::geom_line(
          data = preds,
          ggplot2::aes(x = .data$.x, y = .data$.y, color = .data$group, group = .data$.grp),
          linewidth = 1
        ) +
        ggplot2::labs(x = "Months since birth",
                      y = "Repertoire stability (Jaccard)",
                      color = "Group") +
        theme_phip()

      # normal palette for lines/points
      p <- .add_phip_scales(p, custom_colors)

      # baseline markers
      if (isTRUE(drop_self) && is.data.frame(marker_tbl) && nrow(marker_tbl)) {
        if (!is.null(marker_fill_colors)) {
          # normalize mapping: allow single color, unnamed vector, or named vector
          if (length(marker_fill_colors) == 1L && is.null(names(marker_fill_colors))) {
            marker_fill_colors <- setNames(rep(marker_fill_colors, length(grp_levels)), grp_levels)
          } else if (is.null(names(marker_fill_colors)) &&
                     length(marker_fill_colors) == length(grp_levels)) {
            names(marker_fill_colors) <- grp_levels
          }
          # new fill scale ONLY for markers
          if (!requireNamespace("ggnewscale", quietly = TRUE)) {
            .ph_abort("ggnewscale is required when marker_fill_colors is used.", step = "markers")
          }
          p <- p +
            ggnewscale::new_scale_fill() +
            ggplot2::geom_point(
              data = marker_tbl,
              ggplot2::aes(x = .data$x_marker, y = .data$y_marker, fill = .data$group),
              shape = 21, size = 2.5, stroke = 0.6, colour = "black",
              inherit.aes = FALSE, show.legend = FALSE
            ) +
            ggplot2::scale_fill_manual(values = marker_fill_colors, guide = "none")
        } else {
          # default: use same group palette; hide fill legend
          p <- p +
            ggplot2::geom_point(
              data = marker_tbl,
              ggplot2::aes(x = .data$x_marker, y = .data$y_marker, fill = .data$group),
              shape = 21, size = 2.0, stroke = 0.6, colour = "black",
              inherit.aes = FALSE, show.legend = FALSE
            ) +
            ggplot2::guides(fill = "none")
        }
      }

      p
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ==============================================================================
# Correlate stability trajectories (on <phip_repertoire_stability>) and plot
# ==============================================================================

#' @title Correlation heatmap of repertoire *stability* trajectories (subjects × subjects)
#' @description
#' Given a tibble produced by [compute_repertoire_stability()], compute
#' pairwise correlations between **subjects' stability time series** (0..1)
#' within an optional group and plot a heatmap.
#'
#' Each subject is a vector over time: `similarity(time)` vs their baseline.
#' We correlate those vectors between subjects using Pearson/Spearman with
#' pairwise-complete handling of missing time points.
#'
#' @param stab      A tibble of class `"phip_repertoire_stability"` with cols:
#'                  `subject_id`, `group`, `time`, `similarity`. Baseline rows
#'                  may have been dropped depending on compute settings.
#' @param group     Character scalar or NULL. If provided, filter to this group
#'                  level (matching `stab$group`). If NULL, all groups pooled.
#' @param min_times Integer; require at least this many time points per subject
#'                  to be included (default 2).
#' @param method    "pearson" (default) or "spearman".
#' @param require_shared_times Logical; if TRUE, keep only **timepoints shared by
#'                  at least 2 subjects** (reduces spurious ±1 from singletons).
#' @param sort_by_status Logical; if TRUE, order subjects by status columns.
#' @param subject_meta Optional tibble with `subject_id` plus optional
#'                  `pre_status_col`, `post_status_col` used when sorting.
#' @param pre_status_col,post_status_col Character names in `subject_meta`
#'                  for ordering (defaults: "status_t1","status_t2").
#' @param palette   Diverging palette for fill scale (default "RdYlBu").
#' @param limits    Numeric length-2 for color scale; default c(-1, 1).
#' @return A ggplot object; invisibly a list(list(plot=..., corr_matrix=...,
#'         subjects=..., order_info=...)).
#' @export
ph_plot_stability_corr_heatmap <- function(
    stab,
    group     = NULL,
    min_times = 2,
    method    = c("pearson","spearman"),
    require_shared_times = TRUE,
    sort_by_status = FALSE,
    subject_meta = NULL,
    pre_status_col  = "status_t1",
    post_status_col = "status_t2",
    palette  = "RdYlBu",
    limits   = c(-1, 1)
) {
  method <- match.arg(method)
  .data <- rlang::.data

  if (!inherits(stab, "phip_repertoire_stability")) {
    .ph_abort("Input must be a <phip_repertoire_stability> tibble.",
              step = "class check",
              bullets = "Call compute_repertoire_stability() first.")
  }

  .ph_with_timing(
    headline = "Stability trajectory correlation (subjects × subjects)",
    step     = sprintf("method=%s, min_times=%d", method, min_times),
    expr = {
      # --- 1) Filter to group (optional) -------------------------------------
      stab0 <- stab
      if (!is.null(group)) {
        stab0 <- dplyr::filter(stab0, .data$group == !!group)
        .ph_log_info("Filtering to group", bullets = paste("group =", group))
      }

      # --- 2) Keep subjects with enough time points ---------------------------
      stab1 <- stab0 |>
        dplyr::group_by(.data$subject_id) |>
        dplyr::mutate(.n_t = dplyr::n_distinct(.data$time)) |>
        dplyr::ungroup() |>
        dplyr::filter(.data$.n_t >= !!min_times) |>
        dplyr::select(-.data$.n_t)

      n_keep <- dplyr::n_distinct(stab1$subject_id)
      .chk_cond(n_keep < 2L, "Fewer than 2 subjects remain after filtering; cannot compute correlations.")

      # --- 3) Optionally restrict to shared time points -----------------------
      if (isTRUE(require_shared_times)) {
        shared_times <- stab1 |>
          dplyr::distinct(.data$subject_id, .data$time) |>
          dplyr::count(.data$time, name = "n_subj") |>
          dplyr::filter(.data$n_subj >= 2L) |>
          dplyr::pull(.data$time)
        stab1 <- dplyr::filter(stab1, .data$time %in% shared_times)
        .ph_log_info("Keeping timepoints shared by ≥2 subjects",
                     bullets = sprintf("timepoints kept: %d", length(shared_times)))
      }

      # --- 4) Build subjects × time wide matrix of stability ------------------
      # Rows: subjects; Cols: times; Values: similarity
      wide <- stab1 |>
        dplyr::select(.data$subject_id, .data$time, .data$similarity) |>
        tidyr::pivot_wider(names_from = .data$time,
                           values_from = .data$similarity) |>
        dplyr::arrange(.data$subject_id)

      subj_ids <- wide$subject_id
      mat <- as.matrix(wide[ , setdiff(colnames(wide), "subject_id"), drop = FALSE])

      # --- 5) Correlation across subjects (pairwise complete) -----------------
      # cor between rows -> cor(t(mat))
      if (ncol(mat) == 0L) {
        .ph_abort("No time columns available after filtering.",
                  step = "build matrix")
      }
      C <- try(stats::cor(t(mat), use = "pairwise.complete.obs", method = method), silent = TRUE)
      if (inherits(C, "try-error")) {
        .ph_abort("Correlation failed", step = "stats::cor", bullets = "Check for constant or empty vectors.")
      }
      dimnames(C) <- list(subj_ids, subj_ids)

      # --- 6) Optional ordering by status (metadata) --------------------------
      order_info <- NULL
      if (isTRUE(sort_by_status)) {
        if (is.null(subject_meta)) {
          .ph_warn("sort_by_status=TRUE but subject_meta is NULL; using alphabetical order.")
        } else {
          need_cols <- c("subject_id", pre_status_col, post_status_col)
          miss <- setdiff(need_cols, colnames(subject_meta))
          if (length(miss)) {
            .ph_warn("subject_meta missing columns; falling back to alphabetical.",
                     bullets = paste("missing:", paste(miss, collapse = ", ")))
          } else {
            order_info <- subject_meta |>
              dplyr::filter(.data$subject_id %in% subj_ids) |>
              dplyr::arrange(.data[[pre_status_col]], .data[[post_status_col]], .data$subject_id)
            ord <- rev(order_info$subject_id)  # mimic your colleague’s “slice(rev(row_number()))”
            C <- C[ord, ord, drop = FALSE]
          }
        }
      }

      # --- 7) Melt for plotting ----------------------------------------------
      cor_df <- reshape2::melt(C, varnames = c("Var1","Var2"), value.name = "value")
      cor_df$Var1 <- factor(cor_df$Var1, levels = rownames(C))
      cor_df$Var2 <- factor(cor_df$Var2, levels = colnames(C))

      # --- 8) Plot ------------------------------------------------------------
      p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile(colour = "gray60", linewidth = 0.25) +
        ggplot2::scale_fill_distiller(
          palette  = palette,
          direction = -1,
          limits   = limits,
          na.value = "gray90",
          name     = paste0(stringr::str_to_title(method), "\nCorrelation")
        ) +
        ggplot2::coord_fixed() +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_text(size = 6.5, angle = 45, vjust = 1, hjust = 0.5),
          axis.text.y  = ggplot2::element_text(size = 6.5),
          panel.grid   = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank()
        ) +
        ggplot2::labs(
          x = "Subjects (stability trajectories)",
          y = "Subjects (stability trajectories)"
        )

      invisible(list(
        plot        = p,
        corr_matrix = C,
        subjects    = subj_ids,
        order_info  = order_info
      ))

    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ==============================================================================
# phiper logging utilities (ASCII only; based on the chk and cli packages)
# ==============================================================================
# ---- user-tweakable globals (set via options() in .onLoad or zzz.R) ----------
# options(
#   phiper.log.verbose   = TRUE,
#   phiper.log.time_fmt  = "%Y-%m-%d %H:%M:%S",
#   phiper.log.width     = getOption("width", 80)
# )

.ph_opt <- function(key,
                    default = NULL) {
  getOption(paste0("phiper.log.", key), default)
}

.ph_now <- function() {
  format(Sys.time(), .ph_opt("time_fmt", "%H:%M:%S"))
}

.ph_base_prefix <- function(level = "INFO") {
  sprintf("[%s] %-5s ", .ph_now(), toupper(level)[1])
}

# wraps the text nicely, regardless of the console width
.ph_wrap <- function(text,
                     prefix) {
  w <- .ph_opt("width", getOption("width", 80))
  # strwrap: 'initial' for first line, 'prefix' for continuations
  strwrap(text, width = w, initial = prefix, prefix = strrep(
    " ",
    nchar(prefix)
  ))
}

# Compose multi-depth message lines
# currently the maximal supported log depth is 3:
# headline (required), step (optional), bullets (optional chr vec)
.ph_compose_lines <- function(level,
                              headline,
                              step = NULL,
                              bullets = NULL) {
  base <- .ph_base_prefix(level)
  stepP <- paste0(strrep(" ", nchar(base)), "-> ")
  bullP <- paste0(strrep(" ", nchar(base)), "  - ")

  out <- character(0)
  if (!is.null(headline) && nzchar(headline)) {
    out <- c(out, .ph_wrap(headline, base))
  }
  if (!is.null(step) && nzchar(step)) {
    out <- c(out, .ph_wrap(step, stepP))
  }
  if (length(bullets)) {
    for (b in bullets) {
      if (isTRUE(is.na(b)) || !nzchar(b)) next
      out <- c(out, .ph_wrap(b, bullP))
    }
  }
  out
}

# ---- Public logging helpers --------------------------------------------------
## monitor task progress
.ph_log_info <- function(headline,
                         step = NULL,
                         bullets = NULL,
                         verbose = .ph_opt("verbose", TRUE)) {
  if (!isTRUE(verbose)) {
    return(invisible(character()))
  }
  lines <- .ph_compose_lines("INFO", headline, step, bullets)
  cat(paste0(lines, collapse = "\n"), "\n", sep = "")
  invisible(lines)
}

## monitor task progression
.ph_log_ok <- function(headline,
                       step = NULL,
                       bullets = NULL,
                       verbose = .ph_opt("verbose", TRUE)) {
  if (!isTRUE(verbose)) {
    return(invisible(character()))
  }
  lines <- .ph_compose_lines("OK", headline, step, bullets)
  cat(paste0(lines, collapse = "\n"), "\n", sep = "")
  invisible(lines)
}

# Warnings/errors via chk, but formatted to match the style of the logger
.ph_warn <- function(headline, step = NULL, bullets = NULL, ...) {
  lines <- .ph_compose_lines("WARN", headline, step, bullets)
  msg <- paste(lines, collapse = "\n")
  if (requireNamespace("chk", quietly = TRUE)) {
    chk::wrn(msg, ...) # single string -> respects \n
  } else {
    warning(msg, call. = FALSE, ...) # fallback if chk not installed
  }
  invisible(lines)
}

.ph_abort <- function(headline, step = NULL, bullets = NULL, ...) {
  lines <- .ph_compose_lines("ERROR", headline, step, bullets)
  msg <- paste(lines, collapse = "\n")
  if (requireNamespace("chk", quietly = TRUE)) {
    chk::abort_chk(msg, ...) # single string -> respects \n
  } else {
    stop(msg, call. = FALSE, ...) # fallback if chk not installed
  }
}

# ---- Back-compaibility (vecmatch) --------------------------------------------
# drop-in replacement for older .msg() funct
.msg <- function(verbose, ...) {
  if (isTRUE(verbose)) .ph_log_info(headline = paste(...))
  invisible(NULL)
}

# original conditional helper to not break down older code
# upgraded to the unified phiper style
.chk_cond <- function(condition,
                      error_message,
                      error = TRUE,
                      step = NULL,
                      bullets = NULL,
                      ...) {
  # log nopthing
  if (!isTRUE(condition)) {
    return(invisible(FALSE))
  }

  # print error and abort exec
  if (isTRUE(error)) {
    .ph_abort(headline = error_message, step = step, bullets = bullets, ...)
  } else {
    # print warning and go on
    .ph_warn(headline = error_message, step = step, bullets = bullets, ...)
  }
  invisible(TRUE)
}

# ---- timing helper for sections ----------------------------------------------
## many tasks in phiper can be long/take a while; it was important to have the
## infos on timing - this func wraps a task to get a start/end pair in the same
## style
.ph_with_timing <- function(headline,
                            step = NULL,
                            bullets = NULL,
                            expr,
                            verbose = .ph_opt("verbose", TRUE)) {
  t0 <- Sys.time()
  .ph_log_info(headline = headline, step = step, bullets = bullets, verbose = verbose)

  res <- tryCatch(
    {
      force(expr)
    }, # evaluate user's code
    finally = {
      dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3)
      .ph_log_ok(
        headline = paste0(headline, " - done"),
        step     = sprintf("elapsed: %ss", dt),
        verbose  = verbose
      )
    }
  )

  res
}

# ==============================================================================
# phiper checks + additional helpers (ASCII-only, unified with phiper logger)
# it depends on: .ph_abort(), .ph_warn(), .chk_cond(), word_list(), add_quotes()
# ==============================================================================

# -- check if filename has given extension ------------------------------------
# comes in handy when loading .csv or .parquet; provide filename and vector of
# extensions to check (eg c(".csv", ".parquet"))
.chk_extension <- function(name,
                           x_name,
                           ext_vec) {
  if (is.null(ext_vec) || !length(ext_vec)) {
    return(invisible(TRUE))
  }

  base <- basename(name %||% "") # extracting filename from paths
  parts <- strsplit(base, "\\.", fixed = FALSE)[[1]] # names + complex ext

  # taking last ext after . (eg .tar.gz --> .gz)
  ext <- if (length(parts) > 1L) {
    tolower(paste(parts[-1L], collapse = "."))
  } else {
    ""
  }

  norm <- function(x) sub("^\\.+", "", tolower(x)) # normalize ext
  got <- if (nzchar(ext)) norm(ext) else "<none>"
  ok <- nzchar(ext) && got %in% norm(ext_vec) # final ext check

  if (!ok) {
    .ph_abort(
      headline = sprintf("Invalid file extension for `%s`.", x_name),
      step = sprintf("validating path: %s", name),
      bullets = c(
        sprintf("got: %s", add_quotes(got, 2L)),
        sprintf(
          "allowed: %s",
          word_list(add_quotes(norm(ext_vec), 2L), and_or = "or")
        )
      )
    )
  }
  invisible(TRUE)
}

# -- check if NULL and replace with default when TRUE (warn in unified style) --
.chk_null_default <- function(x,
                              x_name,
                              method,
                              default) {
  if (is.null(x)) {
    # format the default for print
    fmt <- function(v) {
      if (is.character(v) && length(v) == 1L) {
        return(add_quotes(v, 2L))
      }
      if (is.atomic(v) && length(v) == 1L) {
        return(as.character(v))
      }
      sprintf("<%s>", paste(class(v), collapse = "/"))
    }

    # generate warning and the replace
    .ph_warn(
      headline = sprintf("Argument `%s` missing; using default.", x_name),
      step     = sprintf("method: %s", add_quotes(method, 2L)),
      bullets  = sprintf("default: %s", fmt(default))
    )
    x <- default
  }
  x
}

# -- validate path to file -----------------------------------------------------
.chk_path <- function(path,
                      arg_name,
                      extension) {
  ## error when path not a string
  .chk_cond(
    !chk::vld_string(path),
    sprintf("`%s` must be a character scalar.", arg_name),
    step    = "path validation",
    bullets = sprintf("got class: %s", paste(class(path), collapse = "/"))
  )

  ## error when path does not exist
  .chk_cond(
    !chk::vld_file(path),
    sprintf("File for `%s` does not exist.", arg_name),
    step    = "path validation",
    bullets = sprintf("path: %s", path)
  )

  # optionally extension check if provided
  if (!missing(extension) && length(extension)) {
    .chk_extension(
      path,
      arg_name,
      extension
    )
  }

  invisible(TRUE)
}

# -- clean wordlists for message generation ------------------------------------
# for multiple arguments/values
word_list <- function(word_list = NULL,
                      and_or = "and",
                      is_are = FALSE,
                      quotes = FALSE) {
  # Make "a and b" / "a, b, and c"; optionally append "is/are".
  word_list <- setdiff(word_list, c(NA_character_, ""))

  if (is.null(word_list)) {
    out <- ""
    attr(out, "plural") <- FALSE
    return(out)
  }

  word_list <- add_quotes(word_list, quotes)

  len_wl <- length(word_list)

  if (len_wl == 1L) {
    out <- word_list
    if (is_are) out <- paste(out, "is")
    attr(out, "plural") <- FALSE
    return(out)
  }

  if (is.null(and_or) || isFALSE(and_or)) {
    out <- paste(word_list, collapse = ", ")
  } else {
    and_or <- match.arg(and_or, c("and", "or"))
    if (len_wl == 2L) {
      out <- sprintf("%s %s %s", word_list[1L], and_or, word_list[2L])
    } else {
      out <- sprintf(
        "%s, %s %s",
        paste(word_list[-len_wl], collapse = ", "),
        and_or, word_list[len_wl]
      )
    }
  }

  if (is_are) out <- sprintf("%s are", out)
  attr(out, "plural") <- TRUE
  out
}

# -- quoting helper (unified error style) --------------------------------------
# define number of quotes you want --> for printing logs/messages/warnings
# or define the quotes itself as a string
add_quotes <- function(x,
                       quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }
  if (isTRUE(quotes)) quotes <- '"'

  if (chk::vld_string(quotes)) {
    return(paste0(quotes, x, quotes))
  }

  if (!chk::vld_count(quotes) || quotes > 2) {
    .ph_abort(
      headline = "Invalid `quotes` argument.",
      step = "formatting add_quotes()",
      bullets = c(
        "allowed: FALSE, TRUE, 0, 1, 2, or a single-character string",
        sprintf("got class: %s", paste(class(quotes), collapse = "/"))
      )
    )
  }

  if (quotes == 0L) {
    return(x)
  }
  if (quotes == 1L) {
    return(sprintf("'%s'", x))
  }
  sprintf('"%s"', x)
}

# -- not-in operator -----------------------------------------------------------
`%nin%` <- function(x, inx) {
  !(x %in% inx)
}

# -- NULL-coalescing helper ----------------------------------------------------
`%||%` <- function(x, y) if (!is.null(x)) x else y


# ==============================================================================
# database helpers
# ==============================================================================
# ensure peptide_library is queryable from the SAME connection as data_long
.ensure_peplib_on_main <- function(x, schema_alias = "peplib") {
  main_con <- dbplyr::remote_con(x$data_long)
  pep_con <- if (!is.null(x$meta$peptide_con)) {
    x$meta$peptide_con
  } else {
    dbplyr::remote_con(x$peptide_library)
  }

  # try zero-copy ATTACH when both are DuckDB
  if (inherits(main_con, "duckdb_connection") &&
      inherits(pep_con, "duckdb_connection")) {
    pep_db_path <- try(pep_con@driver@dbdir, silent = TRUE)
    if (!inherits(pep_db_path, "try-error") &&
        is.character(pep_db_path) &&
        nzchar(pep_db_path)) {
      # ATTACH (ignore "already attached")
      try(
        DBI::dbExecute(
          main_con,
          sprintf(
            "ATTACH '%s' AS %s;",
            pep_db_path,
            schema_alias
          )
        ),
        silent = TRUE
      )

      # Detect base table name (fallback: "peptide_meta")
      base_name <- tryCatch(
        {
          nm <- dbplyr::remote_name(x$peptide_library)
          if (is.null(nm) || !nzchar(nm)) {
            "peptide_meta"
          } else {
            sub("^.*\\.", "", nm)
          }
        },
        error = function(e) "peptide_meta"
      )

      # Try both two-part and three-part references
      try_tbl <- function(sql_expr) {
        tryCatch(dplyr::tbl(main_con, dbplyr::sql(sql_expr)),
                 error = function(e) NULL
        )
      }

      candidates <- c(
        sprintf("SELECT * FROM %s.%s", schema_alias, base_name),
        sprintf("SELECT * FROM %s.main.%s", schema_alias, base_name)
      )

      for (q in candidates) {
        out <- try_tbl(q)
        if (!is.null(out)) {
          return(out)
        }
      }
      # If both fail, we’ll fall through to the copy_to() fallback below.
    }
  }

  # Fallback: copy peptidelib into main_con as a TEMP table
  peplib_local <- dplyr::collect(x$peptide_library)
  tmp_name <- paste0("peptide_meta_tmp_", as.integer(Sys.time()))
  dplyr::copy_to(main_con, peplib_local, tmp_name,
                 temporary = TRUE,
                 overwrite = TRUE
  )
}
