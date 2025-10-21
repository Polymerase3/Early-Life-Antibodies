## ----------------------------READING PACKAGES---------------------------------
# setting random generator seed
set.seed(16748991)
seed_before <- .Random.seed

# creating vector of necessary packages
packages <- c(
  "stringr",
  "magrittr",
  "data.table",
  "MASS",
  "ggplot2",
  "multilevelTools",
  "dplyr",
  "ggpubr",
  "tidyr",
  "lmerTest",
  "lme4",
  "extraoperators",
  "JWileymisc",
  "showtext",
  "DBI",
  "duckdb",
  "knitr",
  "kableExtra",
  "htmltools",
  "webshot2",
  "magick",
  "ggExtra",
  "mclust",
  "dplyr",
  "DBI",
  "readxl",
  "stringi",
  "glue",
  "purrr",
  "forcats",
  "DescTools",
  "corrplot",
  "tibble",
  "wesanderson",
  "tidytext",
  "Rtsne",
  "Matrix",
  "plotly",
  "randomcoloR",
  "stringr",
  "scales",
  "readr",
  "ggsignif",
  "ggnewscale",
  "scales",
  "ggbreak",
  "RColorBrewer",
  "rgl",
  "htmlwidgets",
  "forcats",
  "writexl"
)

## load development version of phiper
devtools::install_github("Polymerase3/phiper", force = TRUE)
devtools::load_all("/home/noxia/Documents/R/phiper")
Sys.sleep(3)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  pak::pkg_install(packages[!installed_packages])
}

# packages loading
invisible(lapply(packages, library, character.only = TRUE))

# add font + phiper use font in all plots
font_add_google("Montserrat", "monte")
library(phiper)
phip_use_montserrat()
showtext_auto()


# removing unnecessary variables
rm(list = c("installed_packages", "packages"))

## ------------------ PHIP DATA + ages -----------------------------------------
## create phip_data
withr::with_preserve_seed({
  ps <- phip_convert(
    data_long_path = "data/babies_counts.parquet",
    backend = "duckdb",
    peptide_library = TRUE,
    subject_id = "subjectID",
    peptide_id = "peptideID",
    materialise_table = TRUE,
    auto_expand = FALSE,
    n_cores = 10
  )
})

## load the ages data.frame
ages <- read_excel("data/time_collection_blood_LLNEXT.xlsx",
  col_types = c(
    "text", "text", "text",
    "text", "text", "numeric", "text",
    "text", "text", "text", "text", "numeric",
    "text", "numeric", "numeric", "numeric"
  )
)

ages$next_id <- gsub("LLNEXT", "", ages$next_id, fixed = TRUE)

## Converting TIMEPOINTS to factors
ages$timepoint_factor <- dplyr::recode_factor(ages$timepoint,
  "P12" = "T0",
  "P28" = "T1",
  "B" = "T2",
  "W2" = "T3",
  "M1" = "T4",
  "M2" = "T5",
  "M3" = "T6",
  "M6" = "T7",
  "M12" = "T8"
)

## NOTE: after i performed the standardization, there were some observations
## missing for the following subjects:
# subject_id timepoint big_group
# <chr>          <dbl> <chr>
# 1 010577            38 kid_serum
# 2 008101            86 kid_serum
# 3 008525            38 mom_serum
# 4 206131            38 mom_serum
#
# The first two are obvious missing data in the ages data.frame, they have NA's
# in the place of exact age; for the subject 010577 there was no exact_age at
# birth, but we can approximate it with 0; for the subject 008101 exact age at
# 12 months was missing, so i approximated it with adding ~ 9 months to the M3
# value: 122.3229167 + 9 * 30
# For the other two mums, the exact_ages at birth were also missing, i resolved
# it manually using the last given value.

ages$exact_age[ages$next_id == "010577" & ages$timepoint == "B"] <- 0
ages$exact_age[ages$next_id == "008101" &
  ages$timepoint == "M12"] <- 122.3229167 + 9 * 30
ages$exact_age[ages$next_id == "008525" &
  ages$timepoint == "B"] <- 11148.52 + 22 * 7
ages$exact_age[ages$next_id == "206131" &
  ages$timepoint == "B"] <- 13839.48 + 22 * 7

## ----- standardizing the ages to days_since_birth on a common scale ----------
## calculating the mom age at birth
mom_birth_age <- ages %>%
  filter(type == "M") %>%
  group_by(next_id) %>%
  arrange(
    !(timepoint == "B" & !is.na(exact_age)),
    exact_age
  ) %>%
  slice(1) %>%
  transmute(next_id,
    mom_minimal_age = exact_age,
    timepoint
  ) %>%
  ungroup()

## for some moms there is no data at birth --> i have to extrapolate from the
## timepoints
mom_birth_age <- mom_birth_age %>%
  mutate(
    mom_birth_age = case_when(
      timepoint == "B" ~ mom_minimal_age, # nothin changes
      timepoint == "P12" ~ mom_minimal_age + 26 * 7, # 38 weeks = 266 days
      timepoint == "P28" ~ mom_minimal_age + 10 * 7, # 22 weeks = 154 days
      TRUE ~ NA_real_ # fallback
    )
  ) %>%
  dplyr::select(next_id, mom_birth_age)

## merge the data with the ages
withr::with_preserve_seed({
  ps_merged <- left_join(ps,
    ages[, c("next_id", "timepoint_factor", "exact_age")],
    by = dplyr::join_by(
      subject_id == next_id,
      timepoint_factor == timepoint_factor
    ),
    copy = TRUE
  )
})

## selecting only necessary variables
ps_merged %<>% dplyr::select(
  sample_id, subject_id, big_group,
  timepoint, timepoint_factor, peptide_id, relative,
  dyade, exact_age, fold_change
)

## adding mom_birth_age
withr::with_seed(1653, {
  ps_merged <- left_join(ps_merged, mom_birth_age,
    by = dplyr::join_by(subject_id == next_id),
    copy = TRUE
  )
})

# converting the exact_age to days_since_birth variable
ps_merged %<>%
  mutate(
    months_since_birth = case_when(
      big_group == "kid_serum" ~ exact_age / 30,
      big_group == "mom_serum" ~ (exact_age - mom_birth_age) / 30,
      big_group == "mom_milk" & timepoint == 40 ~ 0.5,
      big_group == "mom_milk" & timepoint == 42 ~ 1,
      big_group == "mom_milk" & timepoint == 46 ~ 2,
      big_group == "mom_milk" & timepoint == 50 ~ 3,
      big_group == "mom_milk" & timepoint == 62 ~ 6,
      big_group == "mom_milk" & timepoint == 86 ~ 12,
      # dont have measurments --> approx
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::select(-mom_birth_age, -timepoint, -exact_age)

## sanity check --> any missing values in the timepoint?
ps_merged %>%
  summarise(
    n_missing = sum(is.na(months_since_birth))
  )

## remove unnecessary
rm(list = c("ps", "ages", "mom_birth_age"))
gc()

## ------------------- APPENDING METADATA --------------------------------------
## we will need them anyway for the GLMMs and tSNEs, so better to append them rn
## the metadata are actually only important for the kids, so i will append them
## only to the kids --> we dont do any grouping analyses for moms/milks -->
## wouldn't make much sense actually
metas <- read_tsv(
  file = "data/Metadata_LLNEXT_including_feedingmode.txt",
  na = c("NA", "", "n/a", "N/A"),
  trim_ws = TRUE,
  locale = locale(encoding = "UTF-8"),
  col_types = cols(
    NG_ID = col_character(),
    next_id_infant = col_character(),
    next_id_mother = col_character(),
    FAMILY.x = col_character(),
    infant_relations = col_character(),
    sibling_number = col_character(),
    twin_pair = col_character(),
    birth_deliverybirthcard_mode_binary = col_character(),
    infant_birthcard_feeding_mode_after_birth = col_character(),
    infant_birthcard_med_antibiotics_T195_J01CA_B = col_character(),
    mother_birthcard_age_at_delivery = col_double(),
    infant_misc_sex = col_character(),
    infant_ffq_ever_never_breastfed = col_character(),
    birth_birthcard_reg_gestational_age_days = col_double(),
    birth_birthcard_reg_gestational_age_weeks = col_double(),
    birth_birthcard_reg_gestational_age_week_string = col_character(),
    birth_birthcard_reg_source = col_character(),
    feedmode_m3 = col_character(),
    feedmode_m3_txt = col_character(),
    baby_birthcard_delivery_mode = col_character(),
    Feeding_mode_Birth = col_character(),
    mother_health_smoked_one_whole_year_p18 = col_character(),
    mother_birthcard_parity = col_double(),
    birth_deliverybirthcard_place_delivery_simple = col_character(),
    FAMILY.y = col_character(),
    Timepoint_categorical = col_character(),
    infant_ffq_feeding_mode_simple = col_character(),
    infant_ffq_feeding_mode_complex = col_character(),
    infant_ffq_breastfeeding_freq_weighted = col_double(),
    infant_ffq_formula_follow_on_freq_weighted = col_double(),
    infant_ffq_stopped_breastfeeding = col_character()
  )
) %>%
  # trim any leftover whitespace in character columns
  mutate(across(where(is.character), ~ trimws(.))) %>%
  # tidy duplicated family columns: prefer FAMILY.x, if missing use FAMILY.y
  mutate(
    FAMILY = coalesce(FAMILY.x, FAMILY.y)
  ) %>%
  select(-FAMILY.x, -FAMILY.y) # drop originals

# select the vars of interest: other are either duplicates or not-necessary or
# doesn't add any value (antibiotics)
metas %<>% select(
  NG_ID,
  next_id_infant, # each infant is unique
  birth_deliverybirthcard_mode_binary,
  infant_birthcard_feeding_mode_after_birth,
  mother_birthcard_age_at_delivery,
  infant_misc_sex,
  infant_ffq_ever_never_breastfed,
  birth_birthcard_reg_gestational_age_days,
  birth_birthcard_reg_gestational_age_weeks,
  feedmode_m3,
  mother_health_smoked_one_whole_year_p18,
  mother_birthcard_parity,
  birth_deliverybirthcard_place_delivery_simple,
  Timepoint_categorical,
  infant_ffq_feeding_mode_simple,
  infant_ffq_feeding_mode_complex
)

## the problem in the new metadata is, that there are multiple rows per subject;
## it is because somebody appended the feeding mode at the end of the old meta
## table with multiple timepoints pro subject. We now have to pivot them to wide
id_col <- "next_id_infant"
time_col <- "Timepoint_categorical"
pivot_cols <- c("infant_ffq_feeding_mode_simple",
                "infant_ffq_feeding_mode_complex")

# collapse non-timepoint columns: take first non-NA per subject
other <- metas %>%
  select(-all_of(c(time_col, pivot_cols))) %>%
  group_by(across(all_of(id_col))) %>%
  summarise(across(everything(), ~ first(na.omit(.x))), .groups = "drop")

# pivot timepoint feeding cols wide
time_wide <- metas %>%
  select(all_of(c(id_col, time_col, pivot_cols))) %>%
  distinct() %>%
  pivot_wider(
    id_cols = all_of(id_col),
    names_from = all_of(time_col),
    values_from = all_of(pivot_cols),
    names_sep = "_"
  )

# join -> one row per subject
metas <- other %>%
  left_join(time_wide, by = id_col) %>%
  rename(
    id_infant = next_id_infant,
    NG_ID = NG_ID,
    delivery_mode = birth_deliverybirthcard_mode_binary,
    feedmode_birth = infant_birthcard_feeding_mode_after_birth,
    mother_delivery_age = mother_birthcard_age_at_delivery,
    infant_sex = infant_misc_sex,
    ever_breastfed = infant_ffq_ever_never_breastfed,
    gestage_birth_days = birth_birthcard_reg_gestational_age_days,
    gestage_birth_weeks = birth_birthcard_reg_gestational_age_weeks,
    feedmode_m3 = feedmode_m3,
    smoking = mother_health_smoked_one_whole_year_p18,
    parity = mother_birthcard_parity,
    delivery_place = birth_deliverybirthcard_place_delivery_simple
  )

# shorten the long ffq column prefixes:
names(metas) <- names(metas) %>%
  gsub("^infant_ffq_feeding_mode_simple_", "ffq_simple_", .) %>%
  gsub("^infant_ffq_feeding_mode_complex_", "ffq_complex_", .)

# delete the LLNEXT prefix from the ID
metas$id_infant <- str_sub(metas$id_infant, 7)

# delete the NG_ID
metas %<>% select(-NG_ID)

# convert the character NAs to real NAs for all variables
metas <- metas %>%
  mutate(across(where(is.character), ~ na_if(.x, "NA")))

# duckdb-optimised merge of the metadata to the main ps_merged object;
# only for the infants --> all other data will have NAs in the place of the
# metadata!!!

# same DuckDB connection as the big phiper object
con <- dbplyr::remote_con(ps_merged$data_long)

# 1) copy metas into duckdb (TEMP) so no auto-copy happens
tmp <- paste0("metas_tmp_", as.integer(Sys.time()))
copy_to(con, metas, name = tmp, temporary = TRUE, overwrite = TRUE)
metas_db <- tbl(con, tmp)

# 2) LEFT JOIN directly on the tbl_lazy then MATERIALIZE that tbl
joined_tbl <- left_join(
  ps_merged$data_long,
  metas_db,
  by = dplyr::join_by(subject_id == id_infant)
)

joined_tbl_mat <- compute(
  joined_tbl,
  name = paste0("ps_merged_long_", as.integer(Sys.time())),
  temporary = TRUE
)

# 3) put the materialized table back into phip_data
ps_merged$data_long <- joined_tbl_mat
ps_merged_meta <- ps_merged

# 4) now its safe to drop the temp metas table
DBI::dbExecute(con, sprintf('DROP TABLE IF EXISTS "%s"', tmp))

# 5) clean other mess
rm(
  "ps_merged", "metas", "con", "joined_tbl", "joined_tbl_mat", "metas_db",
  "tmp"
)
gc()

## ------------------- BOX-COX -------------------------------------------------
## transform the fold_changes by group using the box-cox method
ps_merged_box <- withr::with_seed(
  123,
  transform_fold_boxcox(ps_merged_meta,
    by = "big_group",
    confirm = FALSE
  )
)

## save the meta for kids in tmp folder if not already saved --> important for
## tSNEs!
cache_path <- file.path("data", "meta_kids.rds")

if (file.exists(cache_path)) {
  meta_kids <- readRDS(cache_path)
} else {
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)

  meta_cols <- c(
    "dyade", "delivery_mode", "feedmode_birth", "mother_delivery_age",
    "infant_sex", "ever_breastfed", "gestage_birth_days", "gestage_birth_weeks",
    "feedmode_m3", "smoking", "parity", "delivery_place", "ffq_simple_M1",
    "ffq_simple_M2", "ffq_simple_M9", "ffq_simple_M6", "ffq_simple_M12",
    "ffq_simple_M3", "ffq_simple_W2", "ffq_complex_M1", "ffq_complex_M2",
    "ffq_complex_M9", "ffq_complex_M6", "ffq_complex_M12", "ffq_complex_M3",
    "ffq_complex_W2"
  )

  meta_kids <- ps_merged_box$data_long %>%
    filter(big_group == "kid_serum") %>%
    select(subject_id, sample_id, all_of(meta_cols)) %>%
    group_by(subject_id) %>%
    slice_min(order_by = sample_id, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-sample_id) %>%
    collect()

  saveRDS(meta_kids, cache_path)
}

# meta_kids is now in memory and cached at tmp/meta_kids.rds

# clean mess
rm("ps_merged_meta")
gc()

## ------------------- BINARIZING THE DATA -------------------------------------
# right now the data are usable only for fold_change-specific analyses - it
# consists only of the enriched peptides and it is not possible to directly
# compute prevalences and do the repertoire analyses; we have to expand the data
# to the binary form, which causes the nrow() to explode, as we are storing the
# data in the long format; hence the decision to use duckdb and phiper - the
# long format is necessary for plotting/GLMMs and many other analyses
# expand and tag existence
ps_merged_box_bin <- phip_expand_full_grid(ps_merged_box,
  add_exist = TRUE,
  exist_col = "exist"
)

# then some cosmetics: add new column for the timepoints of breast milk, we
# wanted to recode W2 --> M1 and M2 --> M3; then also recode the dyades,
# originally i hardcoded the same dyades for siblings/twins and mother, so all
# 3 subjects had the same dyade ID; it was a little bit problematic in the
# heatmaps the, as the siblings were plotted against each other, and not the
# same subjects, so i recoded the dyade column to a new one, so we have another
# option to try in the heatmaps

# --- SETTINGS: replace with your pairs (each element length-2) ---
pairs_list <- list(
  c("120677", "150788"),
  c("202903", "202915"),
  c("204122", "205217")
)
subjects <- unlist(pairs_list)
new_tbl <- "ps_merged_box_bin_materialized"
match_col <- "subject_id"

# --- connection + detect source table ---
con <- ps_merged_box_bin$meta$con
orig_tbl <- if (is.character(ps_merged_box_bin$data_long)) {
  ps_merged_box_bin$data_long
} else {
  dbplyr::remote_name(ps_merged_box_bin$data_long)
}
ps_tbl <- tbl(con, orig_tbl)

# --- get current subject->dyade mapping (small table) ---
df_map <- ps_tbl %>%
  select(!!sym(match_col), dyade) %>%
  distinct() %>%
  collect()

df_map <- df_map %>%
  mutate(dyade_num = suppressWarnings(as.integer(trimws(as.character(dyade)))))

# --- compute start = current max numeric dyade (0 if none) ---
start <- if (all(is.na(df_map$dyade_num))) {
  0L
} else {
  max(df_map$dyade_num, na.rm = TRUE)
}

start <- ifelse(is.infinite(start) | is.na(start), 0L, as.integer(start))

# --- allocate new dyade numbers by pair, appended at the end ---
# for pair i -> numbers: start + (2*(i-1)) + 1 and +2
pair_allocs <- lapply(seq_along(pairs_list), function(i) {
  base <- start + 2 * (i - 1)
  tibble(
    subject_id = pairs_list[[i]],
    new_dyade  = as.character(base + seq_len(length(pairs_list[[i]])))
  )
})
mapping <- bind_rows(pair_allocs)

# --- write mapping to DB temporary table ---
DBI::dbWriteTable(con, "tmp_pairs_dyade_map", mapping, temporary = TRUE,
                  overwrite = TRUE)

# --- build lazy recoded table and materialize it ---
# create a new table with timepoint_recoded (same rules) and dyade_recoded
# (coalesce)
ps_tbl %>%
  mutate(
    timepoint_recoded = case_when(
      big_group == "mom_milk" & timepoint_factor == "T3" ~ "T4",
      big_group == "mom_milk" & timepoint_factor == "T5" ~ "T6",
      TRUE ~ timepoint_factor
    )
  ) %>%
  left_join(tbl(con, "tmp_pairs_dyade_map"), by = match_col) %>%
  mutate(dyade_recoded = coalesce(new_dyade, dyade)) %>%
  select(everything(), -new_dyade) %>%
  compute(name = new_tbl, temporary = FALSE)

# --- safe swap: rename original -> backup, new -> original ---
backup_name <- paste0(orig_tbl, "_backup_", format(Sys.time(), "%Y%m%d%H%M%S"))
DBI::dbExecute(con, paste0(
  "ALTER TABLE ", DBI::dbQuoteIdentifier(con, orig_tbl),
  " RENAME TO ", DBI::dbQuoteIdentifier(con, backup_name)
))
DBI::dbExecute(con, paste0(
  "ALTER TABLE ", DBI::dbQuoteIdentifier(con, new_tbl),
  " RENAME TO ", DBI::dbQuoteIdentifier(con, orig_tbl)
))

# --- drop other temp tables (keep only the new data_long) ---
all_tables <- DBI::dbListTables(con)
to_keep <- orig_tbl
to_drop <- setdiff(all_tables, to_keep)
if (length(to_drop) > 0) {
  for (t in to_drop) {
    DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS ",
                               DBI::dbQuoteIdentifier(con, t)))
  }
}

# refresh R reference, ramove objects and double-free memory
ps_merged_box_bin$data_long <- dplyr::tbl(con, orig_tbl)
rm(list = c(
  "con", "mapping", "meta_kids", "other", "ps_merged_box", "ps_tbl",
  "time_wide", "all_tables", "backup_name", "cache_path", "id_col",
  "match_col", "new_tbl", "orig_tbl", "pivot_cols", "t",
  "targets", "time_col", "to_drop", "to_keep", "df_map"
))

gc()
gc()
