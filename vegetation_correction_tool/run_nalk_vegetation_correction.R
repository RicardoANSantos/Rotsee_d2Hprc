# ==============================================================================
# Vegetation-corrected d2H of precipitation from sedimentary n-C29 alkanes
# Monte Carlo propagation + GAM smoothing + Grapher-ready exports
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
#   1. Reads an Excel (.xlsx/.xls) or CSV table of n-alkane data.
#   2. Estimates the grass fraction of the leaf-wax source from the
#      long-chain n-alkane relative-abundance (R.A.) index
#      (Schafer et al., 2016):
#          f_GR = (C31 + C33_used) / (C27 + C31 + C33_used)
#      n-C33 is optional (see `missing_C33_policy`): if the column is absent
#      or a value is missing, the index can fall back to C31 / (C27 + C31).
#   3. Applies a vegetation-weighted apparent fractionation and inverts it to
#      precipitation d2H (Santos et al., 2026):
#          epsilon_mix = f_GR * epsilon_GR + (1 - f_GR) * epsilon_WP
#          d2Hprc      = ((d2H_C29 + 1000) / (1 + epsilon_mix / 1000)) - 1000
#   4. Propagates endmember and analytical uncertainty with a Monte Carlo
#      simulation (default 10,000 draws per sample).
#   5. Fits descriptive GAMs (mgcv, REML) to raw n-C29 d2H, f_GR and d2Hprc,
#      computes the first derivative (rate of change) and reports where the
#      change is statistically significant. Significant parts of each curve
#      are exported as NUMBERS (fitted values), so they can be plotted
#      directly as a highlighted line in Grapher / Excel / any plotting tool.
#   6. Optionally propagates age-model uncertainty into the d2Hprc GAM
#      (age-only sensitivity envelope) when age quantile columns are given.
#   7. Writes everything into a tidy, dataset-specific folder:
#
#      outputs/<dataset>/
#        <dataset>_RESULTS.xlsx         main workbook (README sheet first)
#        <dataset>_METHODS_SUMMARY.txt  methods text, diagnostics, caveats
#        01_tables/                     all tables as CSV
#        02_figures/png/, pdf/          publication-style figures
#        03_diagnostics/                GAM diagnostics, QC, excluded rows,
#                                       residual plots, run manifest,
#                                       R session info
#
# HOW TO USE
#   * Edit section 1 (USER SETTINGS) - normally only `input_file` and the
#     column names/age settings need to change.
#   * Run in RStudio (Source) or from a terminal:
#         Rscript run_nalk_vegetation_correction.R  [input_file]  [sheet]
#   * Relative paths are resolved against the folder containing this script
#     (falling back to the current working directory).
#
# REQUIRED INPUT COLUMNS (names configurable in section 1)
#   age                  sample age (e.g. cal ka BP, cal yr BP, or year CE)
#   n-C29 d2H            per mil VSMOW
#   C27, C31             concentrations (any consistent unit, e.g. ug/g dw)
# OPTIONAL INPUT COLUMNS
#   n-C29 d2H SD         analytical 1 SD (per mil); missing -> treated as 0
#   C33                  concentration; missing handled by missing_C33_policy
#   C29 concentration    carried through only
#   CPI                  flagged when below `cpi_review_threshold`
#   ID / depth columns   carried through to every output (`id_columns`)
#   age quantiles        q025, q16, q84, q975 of the age model (age envelope)
#
# DEPENDENCIES
#   mgcv (required). readxl (for Excel input; a base-R fallback reader is
#   included). writexl (for the .xlsx workbook; CSVs are always written).
#
# LIMITATIONS (details in the README and in the METHODS_SUMMARY output)
#   f_GR assumes n-C31/C33 are mainly grass-derived and n-C27 mainly woody.
#   Some trees (e.g. Fraxinus, Acer) also produce n-C31/C33, so f_GR is a
#   wax-source estimate, not a calibrated estimate of catchment grass cover.
#   Constant epsilon_app endmembers are a first-order approach; test their
#   transferability before applying them outside temperate Central Europe.
#
# REFERENCES
#   Santos, R. N., et al. (2026). Paleoceanography and Paleoclimatology, 41,
#     e2025PA005401. https://doi.org/10.1029/2025PA005401
#   Schafer, I. K., et al. (2016). SOIL, 2, 551-564.
#     https://doi.org/10.5194/soil-2-551-2016
#   Sachse, D., et al. (2012). Annu. Rev. Earth Planet. Sci., 40, 221-249.
#     https://doi.org/10.1146/annurev-earth-042711-105535
#   Wood, S. N. (2017). Generalized Additive Models: An Introduction with R,
#     2nd ed. https://doi.org/10.1201/9781315370279
#   Simpson, G. L. (2018). Modelling palaeoecological time series using GAMs.
#     Front. Ecol. Evol., 6, 149. https://doi.org/10.3389/fevo.2018.00149
# ==============================================================================


# ==============================================================================
# 1. USER SETTINGS
# ==============================================================================

# ---- Input / output -----------------------------------------------------------
input_file  <- "example/example_input.xlsx"  # .xlsx, .xls or .csv
input_sheet <- NULL          # NULL = sheet "n-alkanes" if present, else first
output_root <- "outputs"     # results go to <output_root>/<input file name>/

# ---- Column names in the input table ------------------------------------------
col_age     <- "age_ka_bp"
col_d2H     <- "nalk_C29_d2H_permille"
col_d2H_sd  <- "nalk_C29_d2H_sd"          # optional
col_C27     <- "nalk_C27_conc_ug_g_dw"
col_C29     <- "nalk_C29_conc_ug_g_dw"    # optional, carried through only
col_C31     <- "nalk_C31_conc_ug_g_dw"
col_C33     <- "nalk_C33_conc_ug_g_dw"    # optional
col_CPI     <- "CPI_C23_33"               # optional, used only for a QC flag
id_columns  <- c("Sample_ID", "depth_mid_cm")  # carried through if present

# ---- Age scale ----------------------------------------------------------------
# age_is_BP = TRUE : larger numbers are OLDER (ka BP, yr BP).
# age_is_BP = FALSE: larger numbers are YOUNGER (year CE / AD).
# Rates of change are always reported FORWARD IN TIME, so "increase" means
# the value becomes higher towards the present, whatever the age scale.
age_is_BP        <- TRUE
age_unit         <- "ka"               # used in rate units, e.g. "permille per ka"
age_axis_label   <- "Age (cal ka BP)"
reverse_age_axis <- FALSE              # TRUE plots the oldest age on the left

# Optional age-model quantiles (same unit as col_age). If all four columns are
# present, an age-only uncertainty envelope for the d2Hprc GAM is computed.
# Set to NULL to disable.
age_quantile_columns <- c(
  q025 = "age_ka_bp_q025", q16 = "age_ka_bp_q16",
  q84  = "age_ka_bp_q84",  q975 = "age_ka_bp_q975"
)
n_age_realizations <- 201L

# ---- Vegetation correction ----------------------------------------------------
# "zero"   : a missing C33 value (or a missing C33 column) contributes 0 to the
#            R.A. index, i.e. f_GR = C31 / (C27 + C31) for that sample.
# "exclude": samples without a finite C33 value are excluded.
missing_C33_policy <- "zero"

# Apparent fractionation endmembers (per mil), Santos et al. (2026).
epsilon_WP_mean <- -110; epsilon_WP_sd <- 21   # woody plants
epsilon_GR_mean <- -165; epsilon_GR_sd <- 25   # grasses

n_simulations <- 10000L
random_seed   <- 12345L

# Central estimate reported for d2Hprc and used for the d2Hprc GAM:
# "median" (robust to skewed Monte Carlo output) or "mean" (as in
# Santos et al., 2026).
central_estimate <- "median"

cpi_review_threshold <- 3   # samples with CPI below this are flagged

# ---- GAM settings ---------------------------------------------------------------
gam_basis_k                <- NULL  # NULL = automatic; or an integer >= 3
gam_basis_min              <- 5L
gam_basis_max              <- 60L
gam_observations_per_basis <- 5L
gam_k_check_repetitions    <- 400L
gam_prediction_points      <- 1001L
gam_confidence_level       <- 0.95

# ---- Figures ---------------------------------------------------------------------
figure_dpi <- 600


# ==============================================================================
# 2. SETUP (no edits needed below this line)
# ==============================================================================

missing_C33_policy <- match.arg(missing_C33_policy, c("zero", "exclude"))
central_estimate   <- match.arg(central_estimate, c("median", "mean"))

if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("Package 'mgcv' is required. Install it with install.packages('mgcv').")
}

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  }
  for (frame in rev(sys.frames())) {            # source("script.R")
    if (!is.null(frame$ofile)) return(dirname(normalizePath(frame$ofile, winslash = "/")))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
    path <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(path)) return(dirname(normalizePath(path, winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}
script_dir <- get_script_dir()

command_args <- commandArgs(trailingOnly = TRUE)
if (length(command_args) >= 1 && nzchar(command_args[1])) input_file  <- command_args[1]
if (length(command_args) >= 2 && nzchar(command_args[2])) input_sheet <- command_args[2]

is_absolute_path <- function(path) grepl("^([A-Za-z]:[/\\\\]|/|~)", path)
resolve_path <- function(path, must_exist = FALSE) {
  if (is_absolute_path(path)) return(path.expand(path))
  candidates <- c(file.path(script_dir, path), file.path(getwd(), path))
  if (!must_exist) return(candidates[1])
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) hit[1] else candidates[1]
}

input_path <- resolve_path(input_file, must_exist = TRUE)
if (!file.exists(input_path)) stop("Input file not found: ", input_path)
input_path <- normalizePath(input_path, winslash = "/")
input_ext  <- tolower(tools::file_ext(input_path))

dataset_id <- gsub("[^A-Za-z0-9._-]+", "_", tools::file_path_sans_ext(basename(input_path)))
output_dir <- file.path(resolve_path(output_root), dataset_id)
dirs <- list(
  root  = output_dir,
  table = file.path(output_dir, "01_tables"),
  png   = file.path(output_dir, "02_figures", "png"),
  pdf   = file.path(output_dir, "02_figures", "pdf"),
  diag  = file.path(output_dir, "03_diagnostics")
)
if (dir.exists(output_dir)) {
  message("Output folder exists; files with the same name will be overwritten:\n  ",
    normalizePath(output_dir, winslash = "/"))
}
for (p in dirs) dir.create(p, recursive = TRUE, showWarnings = FALSE)

out_path <- function(dir_key, name, ext) {
  file.path(dirs[[dir_key]], paste0(dataset_id, "_", name, ".", ext))
}

# Files are written to a temporary file first and then copied. If the target is
# open/locked (typical for Excel on Windows), a time-stamped copy is saved
# instead and a warning is given, so a run never fails at the very end.
written_files <- character(0)
safe_write <- function(path, writer) {
  temporary <- tempfile(fileext = paste0(".", tools::file_ext(path)))
  writer(temporary)
  copied <- suppressWarnings(file.copy(temporary, path, overwrite = TRUE))
  if (!isTRUE(copied)) {
    fallback <- file.path(dirname(path), paste0(
      tools::file_path_sans_ext(basename(path)), "_",
      format(Sys.time(), "%Y%m%d_%H%M%S"), ".", tools::file_ext(path)))
    if (!isTRUE(file.copy(temporary, fallback, overwrite = FALSE))) {
      stop("Could not write ", path, " (nor a fallback copy).")
    }
    warning("'", basename(path), "' is open or locked; saved as ", basename(fallback))
    path <- fallback
  }
  unlink(temporary)
  written_files <<- c(written_files, normalizePath(path, winslash = "/"))
  invisible(path)
}
write_csv_safe <- function(x, path) {
  safe_write(path, function(p) utils::write.csv(x, p, row.names = FALSE, na = ""))
}
write_lines_safe <- function(lines, path) {
  safe_write(path, function(p) writeLines(enc2utf8(lines), p, useBytes = TRUE))
}

fmt <- function(x, digits = 3) {
  out <- format(signif(x, digits), trim = TRUE, scientific = FALSE, drop0trailing = TRUE)
  out[!is.finite(x)] <- "NA"
  out
}


# ==============================================================================
# 3. READ INPUT
# ==============================================================================

# ---- Minimal base-R .xlsx reader (used only when readxl is not installed) ----
xml_unescape <- function(x) {
  x <- gsub("&lt;", "<", x, fixed = TRUE); x <- gsub("&gt;", ">", x, fixed = TRUE)
  x <- gsub("&quot;", "\"", x, fixed = TRUE); x <- gsub("&apos;", "'", x, fixed = TRUE)
  gsub("&amp;", "&", x, fixed = TRUE)
}
read_xml_text <- function(path) paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "")
xml_attribute <- function(tag, attribute) {
  m <- regmatches(tag, regexec(paste0("\\b", attribute, "=\"([^\"]*)\""), tag, perl = TRUE))[[1]]
  if (length(m) >= 2) xml_unescape(m[2]) else NA_character_
}
excel_column_number <- function(reference) {
  value <- 0L
  for (ch in strsplit(gsub("[0-9]", "", reference), "")[[1]]) value <- value * 26L + match(toupper(ch), LETTERS)
  value
}
extract_tag_text <- function(xml, tag) {
  hits <- regmatches(xml, gregexpr(paste0("<", tag, "(?:\\s[^>]*)?>(.*?)</", tag, ">"), xml, perl = TRUE))[[1]]
  if (!length(hits)) return(character(0))
  xml_unescape(gsub(paste0("</?", tag, "(?:\\s[^>]*)?>"), "", hits, perl = TRUE))
}
xlsx_sheet_info_base <- function(dir) {
  wb   <- read_xml_text(file.path(dir, "xl", "workbook.xml"))
  rels <- read_xml_text(file.path(dir, "xl", "_rels", "workbook.xml.rels"))
  sheets <- regmatches(wb, gregexpr("<sheet\\b[^>]*/>", wb, perl = TRUE))[[1]]
  rel    <- regmatches(rels, gregexpr("<Relationship\\b[^>]*/>", rels, perl = TRUE))[[1]]
  rel_id <- vapply(rel, xml_attribute, "", attribute = "Id")
  rel_tg <- vapply(rel, xml_attribute, "", attribute = "Target")
  target <- sub("^/", "", rel_tg[match(vapply(sheets, xml_attribute, "", attribute = "r:id"), rel_id)])
  target <- ifelse(grepl("^xl/", target), target, file.path("xl", target))
  data.frame(name = vapply(sheets, xml_attribute, "", attribute = "name"),
    target = target, stringsAsFactors = FALSE)
}
read_excel_base <- function(path, sheet) {
  dir <- tempfile("xlsx_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(path, exdir = dir)
  info <- xlsx_sheet_info_base(dir)
  idx <- match(sheet, info$name)
  if (is.na(idx)) stop("Worksheet '", sheet, "' not found.")
  ss_path <- file.path(dir, "xl", "sharedStrings.xml")
  shared <- if (file.exists(ss_path)) {
    ss <- read_xml_text(ss_path)
    items <- regmatches(ss, gregexpr("<si(?:\\s[^>]*)?>.*?</si>", ss, perl = TRUE))[[1]]
    vapply(items, function(it) paste(extract_tag_text(it, "t"), collapse = ""), "")
  } else character(0)
  xml <- read_xml_text(file.path(dir, info$target[idx]))
  cells <- regmatches(xml, gregexpr("(?s)(?:<c\\b[^>]*?/>|<c\\b[^>]*?>.*?</c>)", xml, perl = TRUE))[[1]]
  if (!length(cells)) stop("Worksheet is empty: ", sheet)
  refs <- vapply(cells, xml_attribute, "", attribute = "r")
  rows <- as.integer(gsub("[^0-9]", "", refs))
  cols <- vapply(refs, excel_column_number, integer(1))
  vals <- vapply(cells, function(cell) {
    type <- xml_attribute(cell, "t")
    if (grepl("/>$", cell) || identical(type, "e")) return(NA_character_)
    if (identical(type, "inlineStr")) {
      t <- extract_tag_text(cell, "t"); return(if (length(t)) paste(t, collapse = "") else NA_character_)
    }
    v <- extract_tag_text(cell, "v")
    if (!length(v)) return(NA_character_)
    if (identical(type, "s")) {
      i <- suppressWarnings(as.integer(v[1])) + 1L
      return(if (is.finite(i) && i <= length(shared)) shared[i] else NA_character_)
    }
    if (identical(type, "b")) return(if (v[1] == "1") "TRUE" else "FALSE")
    v[1]
  }, "")
  m <- matrix(NA_character_, max(rows), max(cols))
  m[cbind(rows, cols)] <- vals
  head <- m[1, ]; keep <- !is.na(head) & nzchar(head)
  body <- m[-1, keep, drop = FALSE]
  body <- body[apply(body, 1, function(r) any(!is.na(r) & nzchar(r))), , drop = FALSE]
  out <- as.data.frame(body, stringsAsFactors = FALSE)
  names(out) <- head[keep]
  out
}
list_excel_sheets_base <- function(path) {
  dir <- tempfile("xlsx_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(path, exdir = dir)
  xlsx_sheet_info_base(dir)$name
}

na_strings <- c("", "NA", "NaN", "NAN", "#DIV/0!", "#N/A", "#VALUE!", "n.d.", "nd", "-")
use_readxl <- requireNamespace("readxl", quietly = TRUE)

if (input_ext == "csv") {
  input_sheet <- NA_character_
  d <- utils::read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = na_strings, fileEncoding = "UTF-8-BOM")
  reader_used <- "utils::read.csv"
} else if (input_ext %in% c("xlsx", "xls")) {
  if (!use_readxl && input_ext == "xls") stop("Reading .xls requires the 'readxl' package.")
  sheets <- if (use_readxl) readxl::excel_sheets(input_path) else list_excel_sheets_base(input_path)
  if (is.null(input_sheet)) input_sheet <- if ("n-alkanes" %in% sheets) "n-alkanes" else sheets[1]
  if (!(input_sheet %in% sheets)) {
    stop("Worksheet '", input_sheet, "' not found. Available: ", paste(sheets, collapse = ", "))
  }
  d <- if (use_readxl) {
    as.data.frame(readxl::read_excel(input_path, sheet = input_sheet, na = na_strings),
      stringsAsFactors = FALSE, check.names = FALSE)
  } else read_excel_base(input_path, input_sheet)
  reader_used <- if (use_readxl) paste0("readxl ", utils::packageVersion("readxl")) else "base-R xlsx reader"
} else {
  stop("Unsupported input type '.", input_ext, "'. Use .xlsx, .xls or .csv.")
}

required <- c(col_age, col_d2H, col_C27, col_C31)
missing_required <- setdiff(required, names(d))
if (length(missing_required)) {
  stop("Input is missing required column(s): ", paste(missing_required, collapse = ", "),
    "\nAvailable columns: ", paste(names(d), collapse = ", "),
    "\nAdjust the column names in section 1 of the script.")
}

as_num <- function(x) suppressWarnings(as.numeric(x))
get_optional <- function(column) if (!is.null(column) && column %in% names(d)) as_num(d[[column]]) else rep(NA_real_, nrow(d))

C33_column_present <- !is.null(col_C33) && col_C33 %in% names(d)
CPI_column_present <- !is.null(col_CPI) && col_CPI %in% names(d)
sd_column_present  <- !is.null(col_d2H_sd) && col_d2H_sd %in% names(d)
id_present <- intersect(id_columns, names(d))

# Standardised working table. Output column names are fixed so that every
# dataset produces the same file layout.
x <- data.frame(source_row = seq_len(nrow(d)) + 1L)
for (id in id_present) x[[id]] <- d[[id]]
x$age                   <- as_num(d[[col_age]])
x$nalk_C29_d2H_permille <- as_num(d[[col_d2H]])
x$nalk_C29_d2H_sd       <- get_optional(col_d2H_sd)
x$nalk_C27_conc         <- as_num(d[[col_C27]])
x$nalk_C29_conc         <- get_optional(col_C29)
x$nalk_C31_conc         <- as_num(d[[col_C31]])
x$nalk_C33_conc         <- get_optional(col_C33)
x$CPI                   <- get_optional(col_CPI)

age_q_present <- !is.null(age_quantile_columns) && all(age_quantile_columns %in% names(d))
if (age_q_present) {
  for (q in names(age_quantile_columns)) x[[paste0("age_", q)]] <- as_num(d[[age_quantile_columns[[q]]]])
}


# ==============================================================================
# 4. RELATIVE-ABUNDANCE INDEX AND QUALITY CONTROL
# ==============================================================================

x$C33_missing      <- !is.finite(x$nalk_C33_conc)
x$C33_used_as_zero <- x$C33_missing & missing_C33_policy == "zero"
x$nalk_C33_conc_used_for_RA <- ifelse(x$C33_used_as_zero, 0, x$nalk_C33_conc)
x$RA_C33_handling <- ifelse(x$C33_missing,
  if (missing_C33_policy == "zero") "missing -> treated as 0" else "missing -> excluded",
  ifelse(x$nalk_C33_conc < 0, "negative -> invalid", "measured"))

x$long_chain_RA_sum <- x$nalk_C27_conc + x$nalk_C31_conc + x$nalk_C33_conc_used_for_RA
x$f_GR <- (x$nalk_C31_conc + x$nalk_C33_conc_used_for_RA) / x$long_chain_RA_sum
x$f_WP <- 1 - x$f_GR

x$valid_age <- is.finite(x$age)
x$valid_C33 <- if (missing_C33_policy == "zero") {
  x$C33_missing | (is.finite(x$nalk_C33_conc) & x$nalk_C33_conc >= 0)
} else is.finite(x$nalk_C33_conc) & x$nalk_C33_conc >= 0
x$valid_RA <- with(x,
  is.finite(nalk_C27_conc) & nalk_C27_conc >= 0 &
  is.finite(nalk_C31_conc) & nalk_C31_conc >= 0 & valid_C33 &
  is.finite(long_chain_RA_sum) & long_chain_RA_sum > 0 &
  is.finite(f_GR) & f_GR >= 0 & f_GR <= 1)
x$valid_d2H <- is.finite(x$nalk_C29_d2H_permille)
x$included  <- x$valid_age & x$valid_RA & x$valid_d2H

x$analytical_sd_used      <- ifelse(is.finite(x$nalk_C29_d2H_sd) & x$nalk_C29_d2H_sd > 0, x$nalk_C29_d2H_sd, 0)
x$analytical_sd_available <- x$analytical_sd_used > 0
x$CPI_flag <- ifelse(!is.finite(x$CPI), NA_character_,
  ifelse(x$CPI < cpi_review_threshold,
    paste0("CPI < ", cpi_review_threshold, ": review source/degradation"), "OK"))

reasons <- cbind(
  ifelse(!x$valid_age, "missing/invalid age", NA),
  ifelse(!x$valid_RA, "missing/invalid C27/C31/C33 or f_GR outside 0-1", NA),
  ifelse(!x$valid_d2H, "missing n-C29 d2H", NA))
x$exclusion_reason <- apply(reasons, 1, function(r) {
  r <- r[!is.na(r)]; if (length(r)) paste(r, collapse = "; ") else "included"
})

x$sample_notes <- apply(cbind(
  ifelse(x$C33_used_as_zero & C33_column_present, "C33 missing, used as 0", NA),
  ifelse(!x$analytical_sd_available, "no analytical SD (MC uses 0)", NA),
  ifelse(!is.na(x$CPI_flag) & x$CPI_flag != "OK", x$CPI_flag, NA)), 1,
  function(r) paste(r[!is.na(r)], collapse = "; "))

C33_in_index <- any(x$included & !x$C33_missing)


# ==============================================================================
# 5. MONTE CARLO PROPAGATION
# ==============================================================================

summarize_draws <- function(v, prefix) {
  q <- stats::quantile(v, c(0.025, 0.16, 0.50, 0.84, 0.975), names = FALSE, type = 8, na.rm = TRUE)
  out <- c(mean(v, na.rm = TRUE), stats::sd(v, na.rm = TRUE), q)
  names(out) <- paste0(prefix, c("_mean", "_sd", "_q025", "_q16", "_median", "_q84", "_q975"))
  out
}

set.seed(random_seed)
valid_index <- which(x$included)
if (length(valid_index) < 5) stop("Fewer than 5 samples have valid age, R.A. concentrations and n-C29 d2H.")

mc <- t(vapply(valid_index, function(i) {
  eps_wp  <- stats::rnorm(n_simulations, epsilon_WP_mean, epsilon_WP_sd)
  eps_gr  <- stats::rnorm(n_simulations, epsilon_GR_mean, epsilon_GR_sd)
  eps_mix <- x$f_GR[i] * eps_gr + (1 - x$f_GR[i]) * eps_wp
  d2h_fix <- rep(x$nalk_C29_d2H_permille[i], n_simulations)
  d2h_tot <- if (x$analytical_sd_used[i] > 0) {
    stats::rnorm(n_simulations, x$nalk_C29_d2H_permille[i], x$analytical_sd_used[i])
  } else d2h_fix
  c(summarize_draws(eps_mix, "epsilon_mix"),
    summarize_draws(((d2h_fix + 1000) / (1 + eps_mix / 1000)) - 1000, "d2Hprc_method"),
    summarize_draws(((d2h_tot + 1000) / (1 + eps_mix / 1000)) - 1000, "d2Hprc_total"))
}, numeric(21)))

results <- cbind(x[valid_index, ], mc)
results <- results[order(results$age, results$source_row), ]
row.names(results) <- NULL
central_column <- paste0("d2Hprc_total_", central_estimate)
results$d2Hprc_central <- results[[central_column]]

age <- results$age
excluded <- x[!x$included, c("source_row", id_present, "age", "nalk_C29_d2H_permille",
  "nalk_C27_conc", "nalk_C31_conc", "nalk_C33_conc", "RA_C33_handling",
  "valid_age", "valid_RA", "valid_d2H", "exclusion_reason")]
names(excluded)[names(excluded) == "age"] <- col_age


# ==============================================================================
# 6. GAMs, RATES OF CHANGE AND SIGNIFICANT INTERVALS
# ==============================================================================

z_value   <- stats::qnorm(1 - (1 - gam_confidence_level) / 2)
ci_label  <- paste0(round(100 * gam_confidence_level), "")
time_sign <- if (age_is_BP) -1 else 1   # converts d/d(age) into d/d(time)
rate_suffix <- paste0("per_", gsub("[^A-Za-z0-9]+", "_", age_unit))

choose_initial_k <- function(n, n_unique) {
  max_allowed <- min(as.integer(gam_basis_max), n - 1L, n_unique - 1L)
  if (max_allowed < 3L) stop("Too few observations or distinct ages for a GAM smooth.")
  if (!is.null(gam_basis_k)) return(min(max(3L, as.integer(gam_basis_k)), max_allowed))
  min(max(as.integer(gam_basis_min), floor(min(n, n_unique) / gam_observations_per_basis)), max_allowed)
}

k_diagnostic <- function(model, seed_offset) {
  set.seed(random_seed + seed_offset)
  chk <- tryCatch(mgcv::k.check(model, n.rep = gam_k_check_repetitions), error = function(e) NULL)
  if (is.null(chk) || !nrow(chk)) return(c(k_prime = NA, edf = NA, k_index = NA, p_value = NA))
  c(k_prime = chk[1, "k'"], edf = chk[1, "edf"], k_index = chk[1, "k-index"], p_value = chk[1, "p-value"])
}

fit_gam <- function(y, name, seed_offset) {
  md <- data.frame(age = age, y = y)
  md <- md[is.finite(md$age) & is.finite(md$y), ]
  n_unique <- length(unique(md$age))
  if (nrow(md) < 5 || n_unique < 5) stop("GAM '", name, "': at least 5 observations at distinct ages required.")
  max_allowed <- min(as.integer(gam_basis_max), nrow(md) - 1L, n_unique - 1L)
  k0 <- choose_initial_k(nrow(md), n_unique)
  fit_k <- function(k) mgcv::gam(y ~ s(age, k = k, bs = "tp"), data = md,
    family = stats::gaussian(), method = "REML")
  k <- k0; fit <- fit_k(k); chk <- k_diagnostic(fit, seed_offset)
  refitted <- FALSE
  edf_fraction <- chk[["edf"]] / chk[["k_prime"]]
  if (isTRUE(chk[["p_value"]] < 0.05 && edf_fraction > 0.80) && k < max_allowed) {
    k <- min(max_allowed, max(k + 5L, as.integer(ceiling(k * 1.5))))
    fit <- fit_k(k); chk <- k_diagnostic(fit, seed_offset + 1L); refitted <- TRUE
  }
  res <- stats::residuals(fit, type = "deviance")
  gaps <- diff(sort(unique(md$age)))
  s_tab <- summary(fit)$s.table
  list(model = fit, name = name, data = md, k_initial = k0, k = k, refitted = refitted,
    diagnostics = data.frame(
      series = name, n = nrow(md), unique_ages = n_unique,
      age_span = diff(range(md$age)), median_age_spacing = stats::median(gaps),
      maximum_age_gap = max(gaps),
      initial_basis_k = k0, basis_k = k, basis_refitted = refitted,
      k_prime = chk[["k_prime"]], edf = s_tab[1, "edf"],
      edf_over_k_prime = s_tab[1, "edf"] / chk[["k_prime"]],
      k_index = chk[["k_index"]], k_check_pvalue = chk[["p_value"]],
      REML_smoothing_parameter = fit$sp[1],
      deviance_explained = summary(fit)$dev.expl,
      smooth_F = s_tab[1, "F"], smooth_pvalue = s_tab[1, "p-value"],
      residual_rmse = sqrt(mean(res^2)),
      residual_lag1_acf = if (length(res) >= 3) stats::acf(res, plot = FALSE, lag.max = 1)$acf[2] else NA,
      stringsAsFactors = FALSE))
}

series <- list(
  d2H_C29 = list(y = results$nalk_C29_d2H_permille, unit = "permille",
                 label = "raw n-C29 d2H"),
  f_GR    = list(y = results$f_GR, unit = "fraction",
                 label = "grass fraction f_GR"),
  d2Hprc  = list(y = results$d2Hprc_central, unit = "permille",
                 label = paste0("corrected d2Hprc (MC ", central_estimate, ")"))
)

age_grid <- seq(min(age), max(age), length.out = gam_prediction_points)

gam_rate <- function(model) {
  h  <- max(diff(range(age_grid)) * 1e-7, sqrt(.Machine$double.eps))
  lo <- pmax(age_grid - h, min(age_grid)); hi <- pmin(age_grid + h, max(age_grid))
  X  <- (stats::predict(model, data.frame(age = hi), type = "lpmatrix") -
         stats::predict(model, data.frame(age = lo), type = "lpmatrix")) / (hi - lo)
  deriv <- as.numeric(X %*% stats::coef(model))
  se    <- sqrt(pmax(rowSums((X %*% stats::vcov(model)) * X), 0))
  rate  <- time_sign * deriv                       # forward in time
  data.frame(rate = rate, lower = rate - z_value * se, upper = rate + z_value * se)
}

for (s in names(series)) {
  g <- fit_gam(series[[s]]$y, s, seed_offset = match(s, names(series)) * 100L)
  p <- stats::predict(g$model, data.frame(age = age_grid), se.fit = TRUE)
  r <- gam_rate(g$model)
  code <- ifelse(r$lower > 0, 1L, ifelse(r$upper < 0, -1L, 0L))
  series[[s]]$gam  <- g
  series[[s]]$line <- data.frame(fit = as.numeric(p$fit),
    lower = as.numeric(p$fit - z_value * p$se.fit),
    upper = as.numeric(p$fit + z_value * p$se.fit))
  series[[s]]$rate <- r
  series[[s]]$code <- code
}

gam_diagnostics <- do.call(rbind, lapply(series, function(s) s$gam$diagnostics))
gam_diagnostics$percent_of_record_significant <- vapply(series, function(s) 100 * mean(s$code != 0), 1)
row.names(gam_diagnostics) <- NULL

# ---- GAM table (the sheet you plot in Grapher) --------------------------------
# For every series:
#   <s>_GAM, _lower_95, _upper_95   fitted curve and pointwise confidence band
#   <s>_GAM_signif                  fitted value where the rate of change is
#                                   significant (either direction), else blank
#   <s>_GAM_signif_increase/decrease  same, split by direction (forward in time)
#   <s>_rate_per_<unit> (+ CI)      first derivative, forward in time
#   <s>_signif_code                 +1 increase, -1 decrease, 0 not significant
gam_table <- data.frame(age = age_grid)
for (s in names(series)) {
  L <- series[[s]]$line; R <- series[[s]]$rate; code <- series[[s]]$code
  gam_table[[paste0(s, "_GAM")]]                  <- L$fit
  gam_table[[paste0(s, "_GAM_lower_", ci_label)]] <- L$lower
  gam_table[[paste0(s, "_GAM_upper_", ci_label)]] <- L$upper
  gam_table[[paste0(s, "_GAM_signif")]]           <- ifelse(code != 0, L$fit, NA)
  gam_table[[paste0(s, "_GAM_signif_increase")]]  <- ifelse(code == 1, L$fit, NA)
  gam_table[[paste0(s, "_GAM_signif_decrease")]]  <- ifelse(code == -1, L$fit, NA)
  gam_table[[paste0(s, "_rate_", rate_suffix)]]   <- R$rate
  gam_table[[paste0(s, "_rate_lower_", ci_label)]] <- R$lower
  gam_table[[paste0(s, "_rate_upper_", ci_label)]] <- R$upper
  gam_table[[paste0(s, "_signif_code")]]          <- code
}
gam_table$Notes <- apply(vapply(names(series), function(s) {
  code <- series[[s]]$code; R <- series[[s]]$rate
  ifelse(code == 0, NA_character_, paste0(
    s, ": significant ", ifelse(code == 1, "increase", "decrease"),
    " (rate ", fmt(R$rate), " ", series[[s]]$unit, " per ", age_unit,
    ", ", ci_label, "% CI ", fmt(R$lower), " to ", fmt(R$upper), ")"))
}, character(length(age_grid))), 1, function(r) paste(r[!is.na(r)], collapse = "; "))
names(gam_table)[1] <- col_age

# ---- Contiguous significant intervals -----------------------------------------
sample_ages <- age
significant_intervals <- do.call(rbind, lapply(names(series), function(s) {
  code <- series[[s]]$code; L <- series[[s]]$line; R <- series[[s]]$rate
  runs <- rle(code); ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1L
  keep <- runs$values != 0
  if (!any(keep)) return(NULL)
  do.call(rbind, lapply(which(keep), function(k) {
    idx <- starts[k]:ends[k]
    # order the interval forward in time
    t_idx <- if (age_is_BP) rev(idx) else idx
    a_from <- age_grid[t_idx[1]]; a_to <- age_grid[t_idx[length(t_idx)]]
    data.frame(
      series = s,
      direction = if (runs$values[k] == 1) "increase" else "decrease",
      start_age = a_from, end_age = a_to,
      duration = abs(a_to - a_from),
      GAM_value_at_start = L$fit[t_idx[1]],
      GAM_value_at_end = L$fit[t_idx[length(t_idx)]],
      net_change = L$fit[t_idx[length(t_idx)]] - L$fit[t_idx[1]],
      mean_rate = mean(R$rate[idx]),
      max_abs_rate = R$rate[idx][which.max(abs(R$rate[idx]))],
      n_samples_in_interval = sum(sample_ages >= min(a_from, a_to) & sample_ages <= max(a_from, a_to)),
      touches_record_edge = min(idx) == 1L || max(idx) == length(age_grid),
      stringsAsFactors = FALSE)
  }))
}))
if (is.null(significant_intervals)) {
  significant_intervals <- data.frame(series = character(0), direction = character(0),
    start_age = numeric(0), end_age = numeric(0), duration = numeric(0),
    GAM_value_at_start = numeric(0), GAM_value_at_end = numeric(0), net_change = numeric(0),
    mean_rate = numeric(0), max_abs_rate = numeric(0), n_samples_in_interval = integer(0),
    touches_record_edge = logical(0))
}
significant_intervals$caution <- with(significant_intervals, ifelse(
  n_samples_in_interval < 3 & touches_record_edge, "few samples; at record edge",
  ifelse(n_samples_in_interval < 3, "fewer than 3 samples in interval",
  ifelse(touches_record_edge, "at record edge (edge effects likely)", ""))))
names(significant_intervals)[names(significant_intervals) == "start_age"] <- paste0("start_", col_age)
names(significant_intervals)[names(significant_intervals) == "end_age"]   <- paste0("end_", col_age)
names(significant_intervals)[names(significant_intervals) == "mean_rate"] <- paste0("mean_rate_", rate_suffix)
names(significant_intervals)[names(significant_intervals) == "max_abs_rate"] <- paste0("peak_rate_", rate_suffix)


# ==============================================================================
# 7. OPTIONAL AGE-MODEL SENSITIVITY ENVELOPE (d2Hprc GAM)
# ==============================================================================
# Each realisation is a monotonic chronology obtained by reading every sample's
# age at the same quantile u of its own age distribution (piecewise-linear
# between q025, q16, median, q84, q975). The d2Hprc GAM, fitted on the median
# chronology, is warped onto each realisation. Only age-model uncertainty is
# represented; the GAM confidence band and Monte Carlo uncertainty are not.

age_envelope <- NULL
age_envelope_note <- if (is.null(age_quantile_columns)) "disabled" else
  paste0("skipped (columns not found: ",
    paste(setdiff(age_quantile_columns, names(d)), collapse = ", "), ")")
if (age_q_present) {
  qmat <- cbind(results$age_q025, results$age_q16, results$age, results$age_q84, results$age_q975)
  if (any(!is.finite(qmat))) {
    age_envelope_note <- "skipped (missing age quantiles for included samples)"
  } else {
    qmat <- t(apply(qmat, 1, sort))
    u_knots <- c(0.025, 0.16, 0.50, 0.84, 0.975)
    u_grid  <- seq(0.025, 0.975, length.out = n_age_realizations)
    draws <- vapply(u_grid, function(u) apply(qmat, 1, function(q)
      stats::approx(u_knots, q, xout = u, ties = "ordered")$y), numeric(nrow(qmat)))
    monotonic <- apply(draws, 2, function(a) all(diff(a) >= 0))
    draws <- draws[, monotonic, drop = FALSE]
    lower <- max(apply(draws, 2, min)); upper <- min(apply(draws, 2, max))
    if (ncol(draws) >= 2 && upper > lower) {
      env_grid <- seq(lower, upper, length.out = 701)
      pred <- vapply(seq_len(ncol(draws)), function(j) {
        median_equiv <- stats::approx(draws[, j], age, xout = env_grid, ties = "ordered")$y
        as.numeric(stats::predict(series$d2Hprc$gam$model, data.frame(age = median_equiv)))
      }, numeric(length(env_grid)))
      qs <- t(apply(pred, 1, stats::quantile, c(0.025, 0.25, 0.5, 0.75, 0.975),
        names = FALSE, type = 8, na.rm = TRUE))
      age_envelope <- data.frame(env_grid, qs)
      names(age_envelope) <- c(col_age, "d2Hprc_GAM_age_q025", "d2Hprc_GAM_age_q25",
        "d2Hprc_GAM_age_q50", "d2Hprc_GAM_age_q75", "d2Hprc_GAM_age_q975")
      age_envelope_note <- paste0("computed from ", ncol(draws), " monotonic chronologies (",
        sum(!monotonic), " with age reversals dropped)")
    } else {
      age_envelope_note <- "skipped (too few monotonic chronologies)"
    }
  }
}


# ==============================================================================
# 8. TABLES
# ==============================================================================

observations <- data.frame(results[, id_present, drop = FALSE],
  age = results$age,
  d2H_C29_permille_VSMOW = results$nalk_C29_d2H_permille,
  d2H_C29_SD_permille = results$nalk_C29_d2H_sd,
  f_GR = results$f_GR,
  C33_treated_as_zero = results$C33_used_as_zero,
  epsilon_mix_mean_permille = results$epsilon_mix_mean,
  d2Hprc_central_permille_VSMOW = results$d2Hprc_central,
  d2Hprc_lower_68_permille_VSMOW = results$d2Hprc_total_q16,
  d2Hprc_upper_68_permille_VSMOW = results$d2Hprc_total_q84,
  d2Hprc_lower_95_permille_VSMOW = results$d2Hprc_total_q025,
  d2Hprc_upper_95_permille_VSMOW = results$d2Hprc_total_q975,
  Notes = results$sample_notes,
  check.names = FALSE, stringsAsFactors = FALSE)
names(observations)[names(observations) == "age"] <- col_age
names(observations)[names(observations) == "d2Hprc_central_permille_VSMOW"] <-
  paste0("d2Hprc_", central_estimate, "_permille_VSMOW")
if (CPI_column_present) observations$CPI <- results$CPI

mc_summary <- results[, c("source_row", id_present, "age", "nalk_C29_d2H_permille", "nalk_C29_d2H_sd",
  "analytical_sd_used", "nalk_C27_conc", "nalk_C29_conc", "nalk_C31_conc", "nalk_C33_conc",
  "nalk_C33_conc_used_for_RA", "RA_C33_handling", "long_chain_RA_sum", "f_GR", "f_WP",
  "CPI", "CPI_flag", colnames(mc))]
names(mc_summary)[names(mc_summary) == "age"] <- col_age

qc_summary <- data.frame(metric = c(
  "Source rows", "Rows with valid age", "Rows with valid R.A. index", "Rows with n-C29 d2H",
  "Rows included", "Rows excluded", "C33 column present",
  "Included rows with C33 treated as 0", "Included rows with analytical SD",
  paste0("Included rows with CPI < ", cpi_review_threshold),
  "Minimum age", "Maximum age", "Minimum f_GR", "Median f_GR", "Maximum f_GR",
  paste0("Minimum d2Hprc ", central_estimate), paste0("Median d2Hprc ", central_estimate),
  paste0("Maximum d2Hprc ", central_estimate),
  "Mean 95% interval width, method only (permille)",
  "Mean 95% interval width, method + analytical (permille)",
  "Age-model envelope"),
  value = c(nrow(x), sum(x$valid_age), sum(x$valid_RA), sum(x$valid_d2H),
    nrow(results), nrow(excluded), C33_column_present,
    sum(results$C33_used_as_zero & C33_column_present), sum(results$analytical_sd_available),
    if (CPI_column_present) sum(results$CPI < cpi_review_threshold, na.rm = TRUE) else "no CPI column",
    fmt(min(age), 6), fmt(max(age), 6),
    fmt(min(results$f_GR)), fmt(stats::median(results$f_GR)), fmt(max(results$f_GR)),
    fmt(min(results$d2Hprc_central), 4), fmt(stats::median(results$d2Hprc_central), 4),
    fmt(max(results$d2Hprc_central), 4),
    fmt(mean(results$d2Hprc_method_q975 - results$d2Hprc_method_q025)),
    fmt(mean(results$d2Hprc_total_q975 - results$d2Hprc_total_q025)),
    age_envelope_note),
  stringsAsFactors = FALSE)

run_info <- data.frame(setting = c(
  "Script", "Input file", "Worksheet", "Input reader", "Output folder", "Run time",
  "R version", "mgcv version", "Random seed", "Monte Carlo draws per sample",
  "epsilon_WP mean, SD (permille)", "epsilon_GR mean, SD (permille)",
  "Missing C33 policy", "C33 used in the R.A. index", "Central estimate",
  "Age column", "Age scale", "Rate units", "GAM k setting", "GAM k range",
  "GAM prediction points", "GAM confidence level", "Significance rule", "Age-model envelope"),
  value = c("run_nalk_vegetation_correction.R", input_path,
    ifelse(is.na(input_sheet), "(csv)", input_sheet), reader_used,
    normalizePath(output_dir, winslash = "/"), format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    R.version.string, as.character(utils::packageVersion("mgcv")), random_seed, n_simulations,
    paste(epsilon_WP_mean, epsilon_WP_sd, sep = ", "), paste(epsilon_GR_mean, epsilon_GR_sd, sep = ", "),
    missing_C33_policy,
    if (C33_in_index) "yes (measured values)" else "no - f_GR = C31/(C27+C31)",
    paste0("Monte Carlo ", central_estimate), col_age,
    if (age_is_BP) "BP (larger = older)" else "CE (larger = younger)",
    paste0("per ", age_unit, ", forward in time"),
    if (is.null(gam_basis_k)) "automatic" else gam_basis_k,
    paste(gam_basis_min, "-", gam_basis_max), gam_prediction_points, gam_confidence_level,
    paste0("pointwise ", ci_label, "% CI of the first derivative excludes zero"),
    age_envelope_note),
  stringsAsFactors = FALSE)

# ---- README sheet: what every column means ------------------------------------
readme <- rbind(
  data.frame(Sheet = "General", Column = "",
    Description = c(
      paste0("Dataset: ", dataset_id, ". Created ", format(Sys.time(), "%Y-%m-%d %H:%M"), "."),
      "All d2H values are per mil VSMOW. Blank cells are intentional (see below).",
      paste0("Observations: one row per included sample. d2Hprc central value = Monte Carlo ",
        central_estimate, "; 68% = q16-q84; 95% = q2.5-q97.5 (endmember + analytical uncertainty)."),
      paste0("GAM: smooth curves on a regular grid of ", gam_prediction_points, " ages spanning the record."),
      paste0("HOW SIGNIFICANCE WORKS: at each grid age the first derivative (rate of change) of the GAM is ",
        "computed with its pointwise ", ci_label, "% confidence interval. If the whole interval is above ",
        "zero the curve is rising significantly; if it is below zero it is falling significantly; if it ",
        "contains zero the change is not significant."),
      paste0("*_GAM_signif columns contain the fitted GAM value ONLY where the change is significant and are ",
        "blank elsewhere. Plot them on top of *_GAM (e.g. a thicker line in Grapher): the visible pieces are ",
        "the significant parts of the curve. *_signif_increase / *_signif_decrease split them by direction."),
      paste0("Direction is always FORWARD IN TIME (towards the present): 'increase' = value becomes higher ",
        "towards younger ages", if (age_is_BP) " (i.e. towards smaller BP ages)." else "."),
      "The Notes column on the GAM sheet lists, for each grid age, which series change significantly, with the rate and its CI.",
      "Significant_intervals: contiguous significant segments per series, with start/end age (forward in time), net change, mean and peak rates, number of samples inside, and cautions.",
      "Significance flags are pointwise and exploratory (not simultaneous / not corrected for multiple testing). Check GAM_diagnostics and residual plots before interpretation.",
      "LIMITATION: f_GR assumes n-C31/n-C33 come mainly from grasses, but some trees (e.g. Fraxinus, Acer; Santos et al., 2026) also produce them, and high n-alkane producers can bias f_GR. f_GR is a wax-source estimate, not a calibrated estimate of catchment grass cover.",
      "LIMITATION: the correction depends on calibration datasets (temperate Central Europe) and a constant epsilon_app per endmember; it may not fully remove vegetation effects. Test transferability before applying it in other regions and vegetation settings. See METHODS_SUMMARY.txt, section 6."),
    stringsAsFactors = FALSE),
  data.frame(Sheet = "Observations",
    Column = c(id_present, col_age, "d2H_C29_permille_VSMOW", "d2H_C29_SD_permille", "f_GR",
      "C33_treated_as_zero", "epsilon_mix_mean_permille", paste0("d2Hprc_", central_estimate, "_permille_VSMOW"),
      "d2Hprc_lower/upper_68", "d2Hprc_lower/upper_95", "Notes", if (CPI_column_present) "CPI"),
    Description = c(rep("carried through from input", length(id_present)), "sample age",
      "measured n-C29 d2H", "analytical 1 SD (blank = not given)",
      if (C33_in_index) "grass fraction (C31+C33)/(C27+C31+C33)" else "grass fraction C31/(C27+C31) (no C33 available)",
      "TRUE when C33 was missing and set to 0 in f_GR",
      "mean vegetation-weighted apparent fractionation",
      "vegetation-corrected precipitation d2H (central estimate)",
      "central 68% Monte Carlo interval (q16, q84)", "central 95% Monte Carlo interval (q2.5, q97.5)",
      "sample-level flags (missing C33, no analytical SD, low CPI)",
      if (CPI_column_present) "carbon preference index from input"),
    stringsAsFactors = FALSE),
  do.call(rbind, lapply(names(series), function(s) data.frame(Sheet = "GAM",
    Column = paste0(s, c("_GAM", paste0("_GAM_lower/upper_", ci_label), "_GAM_signif",
      "_GAM_signif_increase", "_GAM_signif_decrease", paste0("_rate_", rate_suffix),
      paste0("_rate_lower/upper_", ci_label), "_signif_code")),
    Description = c(
      paste0("GAM fit of ", series[[s]]$label),
      paste0("pointwise ", ci_label, "% confidence band of the fitted mean"),
      "GAM value where the change is significant (either direction); blank = not significant",
      "GAM value where significantly increasing (forward in time); blank otherwise",
      "GAM value where significantly decreasing (forward in time); blank otherwise",
      paste0("first derivative, ", series[[s]]$unit, " per ", age_unit, ", forward in time"),
      paste0("pointwise ", ci_label, "% CI of the rate"),
      "+1 significant increase, -1 significant decrease, 0 not significant"),
    stringsAsFactors = FALSE))),
  data.frame(Sheet = c("GAM", "MC_summary", "GAM_diagnostics", "QC_summary", "Excluded_rows",
    if (!is.null(age_envelope)) "Age_envelope", "Run_info"),
    Column = "",
    Description = c("Notes: plain-language list of significant changes at that grid age",
      "full Monte Carlo summaries (mean, SD, q2.5, q16, median, q84, q97.5) for epsilon_mix, method-only and total d2Hprc",
      "basis size, k-check, edf, deviance explained, residual RMSE and lag-1 autocorrelation per GAM",
      "counts and ranges for this run", "input rows not used and why",
      if (!is.null(age_envelope)) "age-model-only envelope of the d2Hprc GAM (quantiles across chronologies)",
      "all settings used"),
    stringsAsFactors = FALSE)
)


# ==============================================================================
# 9. WRITE TABLES AND WORKBOOK
# ==============================================================================

write_csv_safe(observations,          out_path("table", "observations", "csv"))
write_csv_safe(gam_table,             out_path("table", "GAM_curves", "csv"))
write_csv_safe(significant_intervals, out_path("table", "GAM_significant_intervals", "csv"))
write_csv_safe(mc_summary,            out_path("table", "MonteCarlo_summary", "csv"))
if (!is.null(age_envelope)) write_csv_safe(age_envelope, out_path("table", "d2Hprc_GAM_age_envelope", "csv"))

write_csv_safe(gam_diagnostics, out_path("diag", "GAM_diagnostics", "csv"))
write_csv_safe(qc_summary,      out_path("diag", "QC_summary", "csv"))
write_csv_safe(excluded,        out_path("diag", "excluded_rows", "csv"))
write_csv_safe(run_info,        out_path("diag", "run_settings", "csv"))

workbook_sheets <- list(README = readme, Observations = observations, GAM = gam_table,
  Significant_intervals = significant_intervals, MC_summary = mc_summary,
  GAM_diagnostics = gam_diagnostics, QC_summary = qc_summary, Excluded_rows = excluded)
if (!is.null(age_envelope)) workbook_sheets$Age_envelope <- age_envelope
workbook_sheets$Run_info <- run_info

xlsx_path <- file.path(output_dir, paste0(dataset_id, "_RESULTS.xlsx"))
if (requireNamespace("writexl", quietly = TRUE)) {
  safe_write(xlsx_path, function(p) writexl::write_xlsx(workbook_sheets, p))
} else {
  message("Package 'writexl' not installed: CSV files were written, but not the Excel workbook.\n",
    "Install once with install.packages('writexl').")
}


# ==============================================================================
# 10. FIGURES
# ==============================================================================

pal <- list(blue = "#2B6F9F", blue_dark = "#184A66", green = "#4C7A5B", green_dark = "#2F5238",
  ink = "#1C252B", muted = "#5B6770", grey_dark = "#2E3A42",
  grid_major = "#D8DEE3", grid_minor = "#EDF0F2")
ribbon      <- grDevices::adjustcolor(pal$blue, alpha.f = 0.14)
ribbon_age  <- grDevices::adjustcolor(pal$blue, alpha.f = 0.22)
interval_95 <- grDevices::adjustcolor(pal$ink, alpha.f = 0.30)
interval_68 <- grDevices::adjustcolor(pal$ink, alpha.f = 0.65)
permil <- if (isTRUE(l10n_info()$`UTF-8`)) "‰" else "per mil"

padded_range <- function(v, f = 0.04) {
  r <- range(v, finite = TRUE); s <- diff(r)
  if (!is.finite(s) || s == 0) s <- max(abs(r), 1)
  r + c(-1, 1) * s * f
}
axis_breaks <- function(limits, n = 7) {
  limits <- sort(limits)
  major <- pretty(limits, n = n); major <- major[major >= limits[1] & major <= limits[2]]
  if (length(major) < 2) major <- limits
  step <- min(diff(major)) / 2
  minor <- seq(floor(limits[1] / step) * step, ceiling(limits[2] / step) * step, by = step)
  list(major = major, minor = minor[minor >= limits[1] & minor <= limits[2]])
}

x_limits <- padded_range(age, 0.01)
if (reverse_age_axis) x_limits <- rev(x_limits)
x_ticks <- axis_breaks(x_limits, 8)

new_panel <- function(y_limits) {
  y_ticks <- axis_breaks(y_limits, 6)
  graphics::plot(NA, NA, type = "n", xlim = x_limits, ylim = y_limits,
    axes = FALSE, ann = FALSE, xaxs = "i", yaxs = "i")
  graphics::abline(v = x_ticks$minor, col = pal$grid_minor, lwd = 0.45)
  graphics::abline(h = y_ticks$minor, col = pal$grid_minor, lwd = 0.45)
  graphics::abline(v = x_ticks$major, col = pal$grid_major, lwd = 0.65)
  graphics::abline(h = y_ticks$major, col = pal$grid_major, lwd = 0.65)
  y_ticks
}
finish_axes <- function(y_ticks, x_labels) {
  graphics::axis(1, at = x_ticks$major, labels = if (x_labels) format(x_ticks$major, trim = TRUE) else FALSE,
    tck = -0.020, lwd = 0, lwd.ticks = 0.75, cex.axis = 0.82, col.ticks = pal$ink, col.axis = pal$ink)
  graphics::axis(1, at = setdiff(x_ticks$minor, x_ticks$major), labels = FALSE,
    tck = -0.011, lwd = 0, lwd.ticks = 0.60, col.ticks = pal$ink)
  graphics::axis(2, at = y_ticks$major, labels = format(y_ticks$major, trim = TRUE),
    tck = -0.020, lwd = 0, lwd.ticks = 0.75, cex.axis = 0.82, col.ticks = pal$ink, col.axis = pal$ink)
  graphics::axis(2, at = setdiff(y_ticks$minor, y_ticks$major), labels = FALSE,
    tck = -0.011, lwd = 0, lwd.ticks = 0.60, col.ticks = pal$ink)
  graphics::box(col = pal$ink, lwd = 0.80)
}
draw_gam <- function(s, col, col_sig, band = TRUE) {
  L <- series[[s]]$line
  if (band) graphics::polygon(c(age_grid, rev(age_grid)), c(L$lower, rev(L$upper)), col = ribbon, border = NA)
  graphics::lines(age_grid, L$fit, col = col, lwd = 1.2)
  # significant parts of the curve drawn thicker (same values as *_GAM_signif)
  graphics::lines(age_grid, ifelse(series[[s]]$code != 0, L$fit, NA), col = col_sig, lwd = 3.2)
}
panel_label <- function(txt) graphics::mtext(txt, side = 3, line = 0.25, adj = 0, font = 2, cex = 0.95)

draw_raw_panel <- function(label = "A") {
  sdv <- results$analytical_sd_used
  yt <- new_panel(padded_range(c(results$nalk_C29_d2H_permille + c(-1, 1) %o% sdv, series$d2H_C29$line$lower,
    series$d2H_C29$line$upper)))
  graphics::segments(age, results$nalk_C29_d2H_permille - sdv, age, results$nalk_C29_d2H_permille + sdv,
    col = interval_68, lwd = 0.9)
  draw_gam("d2H_C29", pal$muted, pal$grey_dark)
  graphics::points(age, results$nalk_C29_d2H_permille, pch = 21, bg = "white", col = pal$ink, lwd = 0.65, cex = 0.6)
  finish_axes(yt, FALSE)
  graphics::mtext(bquote(delta^2 * H[italic(n) * "-C"[29]] ~ "(" * .(permil) ~ "VSMOW)"),
    side = 2, line = 3.05, cex = 0.82, las = 0)
  panel_label(label)
}
draw_fgr_panel <- function(label = "B") {
  yt <- new_panel(padded_range(c(results$f_GR, series$f_GR$line$lower, series$f_GR$line$upper)))
  draw_gam("f_GR", pal$green, pal$green_dark)
  graphics::points(age, results$f_GR, pch = 21, bg = pal$green, col = "white", lwd = 0.55, cex = 0.6)
  finish_axes(yt, FALSE)
  lab <- if (C33_in_index) {
    expression(f[GR] == frac(italic(n) * "-C"[31] + italic(n) * "-C"[33],
      italic(n) * "-C"[27] + italic(n) * "-C"[31] + italic(n) * "-C"[33]))
  } else {
    expression(f[GR] == frac(italic(n) * "-C"[31], italic(n) * "-C"[27] + italic(n) * "-C"[31]))
  }
  graphics::mtext(lab, side = 2, line = 2.25, cex = 0.68, las = 0)
  panel_label(label)
}
draw_prc_panel <- function(label = "C", with_age_envelope = FALSE) {
  yv <- c(results$d2Hprc_total_q025, results$d2Hprc_total_q975, series$d2Hprc$line$lower, series$d2Hprc$line$upper)
  if (with_age_envelope) yv <- c(yv, age_envelope[[2]], age_envelope[[6]])
  yt <- new_panel(padded_range(yv))
  if (with_age_envelope) {
    ag <- age_envelope[[1]]
    graphics::polygon(c(ag, rev(ag)), c(age_envelope[[2]], rev(age_envelope[[6]])), col = ribbon, border = NA)
    graphics::polygon(c(ag, rev(ag)), c(age_envelope[[3]], rev(age_envelope[[5]])), col = ribbon_age, border = NA)
  }
  graphics::segments(age, results$d2Hprc_total_q025, age, results$d2Hprc_total_q975, col = interval_95, lwd = 0.65)
  graphics::segments(age, results$d2Hprc_total_q16, age, results$d2Hprc_total_q84, col = interval_68, lwd = 1.15)
  draw_gam("d2Hprc", pal$blue_dark, pal$blue_dark, band = !with_age_envelope)
  graphics::points(age, results$d2Hprc_central, pch = 21, bg = pal$blue, col = "white", lwd = 0.55, cex = 0.6)
  finish_axes(yt, TRUE)
  graphics::mtext(age_axis_label, side = 1, line = 2.65, cex = 0.92)
  graphics::mtext(bquote(delta^2 * H[prc] ~ "(" * .(permil) ~ "VSMOW)"), side = 2, line = 3.05, cex = 0.92, las = 0)
  if (nzchar(label)) panel_label(label)
}

base_par <- function(mar) graphics::par(mar = mar, mgp = c(2.8, 0.72, 0), tcl = -0.22, las = 1,
  lend = "round", xaxs = "i", yaxs = "i", family = "sans")

draw_three_panel <- function() {
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::layout(matrix(1:3, ncol = 1), heights = c(1, 1.15, 1.25))
  base_par(c(0.8, 5.2, 1.4, 1.0)); graphics::par(oma = c(3.7, 0, 0.2, 0))
  draw_raw_panel(); draw_fgr_panel(); draw_prc_panel()
}
draw_single_panel <- function() {
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  base_par(c(4.2, 5.0, 1.2, 1.0)); draw_prc_panel(label = "")
}
draw_age_envelope_panel <- function() {
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  base_par(c(4.2, 5.0, 1.2, 1.0)); draw_prc_panel(label = "", with_age_envelope = TRUE)
}
draw_rate_panels <- function() {
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::layout(matrix(1:3, ncol = 1))
  base_par(c(0.8, 5.2, 1.4, 1.0)); graphics::par(oma = c(3.7, 0, 0.2, 0))
  for (s in names(series)) {
    R <- series[[s]]$rate
    yt <- new_panel(padded_range(c(R$lower, R$upper, 0)))
    graphics::polygon(c(age_grid, rev(age_grid)), c(R$lower, rev(R$upper)), col = ribbon, border = NA)
    graphics::abline(h = 0, col = pal$ink, lwd = 0.7, lty = 2)
    graphics::lines(age_grid, R$rate, col = pal$blue_dark, lwd = 1.2)
    graphics::lines(age_grid, ifelse(series[[s]]$code != 0, R$rate, NA), col = pal$blue_dark, lwd = 3.2)
    finish_axes(yt, s == "d2Hprc")
    graphics::mtext(paste0(s, " rate\n(", series[[s]]$unit, " per ", age_unit, ")"),
      side = 2, line = 2.9, cex = 0.72, las = 0)
    panel_label(LETTERS[match(s, names(series))])
  }
  graphics::mtext(paste0(age_axis_label, "   |   rates forward in time; thick = significant"),
    side = 1, line = 2.3, outer = TRUE, cex = 0.8)
}
draw_residual_diagnostics <- function() {
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::par(mfrow = c(3, 4), mar = c(3.6, 3.8, 3.0, 0.8), mgp = c(2.2, 0.6, 0), cex = 0.7, las = 1,
    cex.main = 0.9)
  for (s in names(series)) {
    m <- series[[s]]$gam$model; r <- stats::residuals(m, type = "deviance"); f <- stats::fitted(m)
    a <- series[[s]]$gam$data$age
    graphics::plot(f, r, pch = 16, cex = 0.6, col = pal$muted, xlab = "fitted", ylab = "residual",
      main = paste(s, "- residuals vs fitted")); graphics::abline(h = 0, lty = 2)
    stats::qqnorm(r, pch = 16, cex = 0.6, col = pal$muted, main = paste(s, "- normal QQ")); stats::qqline(r)
    graphics::plot(a, r, pch = 16, cex = 0.6, col = pal$muted, xlab = col_age, ylab = "residual",
      main = paste(s, "- residuals vs age"), xlim = x_limits); graphics::abline(h = 0, lty = 2)
    stats::acf(r, main = paste(s, "- residual ACF"), na.action = stats::na.pass)
  }
}

render <- function(name, drawer, width, height, dir_png = "png", dir_pdf = "pdf") {
  safe_write(out_path(dir_png, name, "png"), function(p) {
    grDevices::png(p, width = width, height = height, units = "in", res = figure_dpi, bg = "white")
    on.exit(grDevices::dev.off()); drawer()
  })
  safe_write(out_path(dir_pdf, name, "pdf"), function(p) {
    grDevices::pdf(p, width = width, height = height, useDingbats = FALSE)
    on.exit(grDevices::dev.off()); drawer()
  })
}

render("three_panel_d2H_fGR_d2Hprc", draw_three_panel, 6.9, 8.6)
render("d2Hprc", draw_single_panel, 6.9, 3.9)
render("GAM_rates_of_change", draw_rate_panels, 6.9, 7.5)
if (!is.null(age_envelope)) render("d2Hprc_age_model_envelope", draw_age_envelope_panel, 6.9, 3.9)
render("GAM_residual_diagnostics", draw_residual_diagnostics, 11, 8.5, dir_png = "diag", dir_pdf = "diag")


# ==============================================================================
# 11. METHODS SUMMARY, MANIFEST AND SESSION INFO
# ==============================================================================

interval_lines <- if (nrow(significant_intervals)) {
  si <- significant_intervals
  paste0("  ", si$series, ": significant ", si$direction, " from ", fmt(si[[paste0("start_", col_age)]], 4),
    " to ", fmt(si[[paste0("end_", col_age)]], 4), " (", col_age, "); net change ", fmt(si$net_change),
    " ", vapply(si$series, function(s) series[[s]]$unit, ""), "; ", si$n_samples_in_interval, " samples",
    ifelse(nzchar(si$caution), paste0(" [", si$caution, "]"), ""))
} else "  No significant changes were detected in any series."

methods_lines <- c(
  "VEGETATION-CORRECTED d2H PRECIPITATION RECONSTRUCTION - METHODS SUMMARY",
  "======================================================================",
  paste0("Dataset: ", dataset_id),
  paste0("Input: ", input_path, if (!is.na(input_sheet)) paste0(" [sheet: ", input_sheet, "]")),
  paste0("Run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " | ", R.version.string,
    " | mgcv ", utils::packageVersion("mgcv")),
  paste0("Included samples: ", nrow(results), " of ", nrow(x),
    " | age range ", fmt(min(age), 5), " to ", fmt(max(age), 5), " (", col_age, ")"),
  "",
  "0. OUTPUT FILES",
  paste0("  ", dataset_id, "_RESULTS.xlsx   main workbook; start with the README sheet"),
  "  01_tables/        CSV copies of Observations, GAM curves, significant intervals,",
  "                    Monte Carlo summary (and age envelope, if computed)",
  paste0("  02_figures/       PNG (", figure_dpi, " dpi) and PDF figures"),
  "  03_diagnostics/   GAM diagnostics, residual plots, QC summary, excluded rows,",
  "                    run settings, run manifest, R session info",
  "",
  "1. VEGETATION CORRECTION",
  "The grass fraction of the wax source is estimated with the long-chain n-alkane",
  "relative-abundance index (Schafer et al., 2016):",
  "  f_GR = (C31 + C33_used) / (C27 + C31 + C33_used);  f_WP = 1 - f_GR",
  paste0("Missing C33 policy: '", missing_C33_policy, "'. C33 column present: ", C33_column_present, "."),
  if (missing_C33_policy == "zero") paste0("Missing C33 is set to 0 for the index (f_GR = C31/(C27+C31) for those ",
    "samples); ", sum(results$C33_used_as_zero), " included samples use this fallback.")
  else "Samples without a finite C33 value are excluded.",
  "Negative concentrations are always invalid. The index is a wax-source proxy and is",
  "not interchangeable with pollen-based vegetation cover without validation.",
  "",
  "2. FRACTIONATION AND MONTE CARLO PROPAGATION",
  paste0("  epsilon_WP = ", epsilon_WP_mean, " +/- ", epsilon_WP_sd, " permille (normal)"),
  paste0("  epsilon_GR = ", epsilon_GR_mean, " +/- ", epsilon_GR_sd, " permille (normal)"),
  "  epsilon_mix = f_GR * epsilon_GR + (1 - f_GR) * epsilon_WP",
  "  d2Hprc = ((d2H_C29 + 1000) / (1 + epsilon_mix / 1000)) - 1000",
  paste0(n_simulations, " draws per sample (seed ", random_seed, "). 'Method' uncertainty samples the"),
  "endmembers only; 'total' also samples measured n-C29 d2H with its analytical SD",
  "(missing SD = 0 for that component; flagged in the Notes column).",
  paste0("Central estimate: Monte Carlo ", central_estimate, "; 68% interval = q16-q84; 95% = q2.5-q97.5."),
  if (central_estimate == "median") "Santos et al. (2026) report the Monte Carlo mean; state this when comparing directly." else "",
  "",
  "3. GENERALIZED ADDITIVE MODELS",
  "Descriptive GAMs are fitted separately to raw n-C29 d2H, f_GR and the d2Hprc central",
  "estimate (mgcv::gam; Gaussian, identity link, thin-plate spline, REML smoothing).",
  paste0("Initial k = floor(min(n, distinct ages) / ", gam_observations_per_basis, "), limited to ",
    gam_basis_min, "-", gam_basis_max, " and below n. mgcv::k.check tests basis adequacy; if"),
  "p < 0.05 and edf > 80% of k', the model is refitted once with a larger k.",
  paste0("Bands are pointwise ", ci_label, "% confidence intervals of the fitted mean (not prediction"),
  "intervals and not the Monte Carlo interval).",
  "",
  "4. SIGNIFICANT CHANGE - HOW TO READ THE NUMBERS",
  "The first derivative of each GAM (rate of change) is computed on the prediction grid",
  paste0("from the model's linear predictor matrix, with a pointwise ", ci_label, "% confidence interval"),
  "(Simpson, 2018). Rates are expressed FORWARD IN TIME (towards the present):",
  "  - CI entirely above 0  -> significant increase (value rises towards the present)",
  "  - CI entirely below 0  -> significant decrease (value falls towards the present)",
  "  - CI includes 0        -> change not significant",
  "In the GAM sheet/CSV, the columns *_GAM_signif, *_GAM_signif_increase and",
  "*_GAM_signif_decrease hold the GAM fitted value ONLY where the change is significant",
  "and are blank elsewhere. Any number in these columns therefore marks a significant part",
  "of the curve; plotting them over *_GAM highlights those segments (blank cells break the",
  "line in Grapher/Excel). The Notes column states the direction, rate and CI per grid age.",
  "These are pointwise, exploratory flags: they are not simultaneous intervals and are not",
  "corrected for multiple testing, and they do not include chronology uncertainty.",
  "",
  "Significant intervals in this run:",
  interval_lines,
  "",
  "GAM diagnostics for this run:",
  utils::capture.output(print(gam_diagnostics[, c("series", "n", "basis_k", "edf", "k_index",
    "k_check_pvalue", "deviance_explained", "smooth_pvalue", "residual_lag1_acf",
    "percent_of_record_significant")], row.names = FALSE, digits = 3)),
  "",
  "Interpretation checks:",
  "  - edf close to k' together with a low k-check p-value suggests k is too small.",
  "  - Strong lag-1 residual autocorrelation violates the independence assumption and makes",
  "    bands and significance flags too optimistic; consider a GAMM with AR(1)/CAR(1) errors.",
  "  - Intervals with few samples or at the record edges (see 'caution') deserve extra care.",
  "  - Do not extrapolate beyond the observed age range.",
  "",
  "5. AGE-MODEL UNCERTAINTY",
  paste0("Age-model envelope: ", age_envelope_note, "."),
  "When computed, every sample's age is read at the same quantile u of its age distribution",
  "(piecewise-linear between q2.5, q16, median, q84, q97.5) to create monotonic chronologies;",
  "the d2Hprc GAM is warped onto each one and the 2.5/25/50/75/97.5% quantiles are reported.",
  "",
  "6. LIMITATIONS AND APPROPRIATE USE",
  "R.A. index: f_GR assumes that n-C31 and n-C33 come mainly from grasses and n-C27 mainly",
  "from woody plants. Some trees (e.g. Fraxinus, Acer; see Santos et al., 2026) also produce",
  "substantial n-C31 and n-C33, so species with high n-alkane production can bias f_GR even",
  "when they are a small part of the vegetation. f_GR describes the wax source; it is NOT a",
  "calibrated estimate of grass cover in the catchment. The index can also be affected by",
  "degradation, transport and aquatic inputs (check CPI, Paq).",
  "Transferability: vegetation corrections depend on their calibration datasets and may not",
  "fully remove vegetation effects. The endmembers are calibrated for temperate Central",
  "Europe; test transferability before applying them to other regions and vegetation settings.",
  "Constant epsilon_app: constant endmember values are a transparent first-order approach",
  "for assessing the direction and approximate magnitude of vegetation effects. epsilon_app",
  "may not be constant; the Monte Carlo spreads account for part of this, but residual",
  "physiological variability remains an additional source of uncertainty.",
  "",
  "7. REFERENCES",
  "Santos, R. N., Nelson, D. B., Klatt, A., Schubert, C. J., Dubois, N., De Jonge, C., &",
  "  Ladd, S. N. (2026). Central European hydroclimate since the Younger Dryas inferred from",
  "  vegetation-corrected sedimentary plant wax d2H values. Paleoceanography and",
  "  Paleoclimatology, 41, e2025PA005401. https://doi.org/10.1029/2025PA005401",
  "Schafer, I. K., et al. (2016). Leaf waxes in litter and topsoils along a European transect.",
  "  SOIL, 2, 551-564. https://doi.org/10.5194/soil-2-551-2016",
  "Sachse, D., et al. (2012). Annual Review of Earth and Planetary Sciences, 40, 221-249.",
  "  https://doi.org/10.1146/annurev-earth-042711-105535",
  "Simpson, G. L. (2018). Modelling palaeoecological time series using generalised additive",
  "  models. Frontiers in Ecology and Evolution, 6, 149. https://doi.org/10.3389/fevo.2018.00149",
  "Wood, S. N. (2017). Generalized Additive Models: An Introduction with R (2nd ed.).",
  "  Chapman and Hall/CRC. https://doi.org/10.1201/9781315370279"
)
write_lines_safe(methods_lines, file.path(output_dir, paste0(dataset_id, "_METHODS_SUMMARY.txt")))
write_lines_safe(utils::capture.output(utils::sessionInfo()), out_path("diag", "R_sessionInfo", "txt"))

manifest_path <- out_path("diag", "run_manifest", "txt")
write_lines_safe(c(
  paste0(run_info$setting, ": ", run_info$value), "", "Files written:",
  paste0("  ", sub(paste0("^", normalizePath(output_dir, winslash = "/"), "/?"), "", written_files),
    "  (", format(file.size(written_files), big.mark = ","), " bytes)")),
  manifest_path)


# ==============================================================================
# 12. CONSOLE SUMMARY
# ==============================================================================

cat("\n== n-alkane vegetation correction ==\n")
cat("Input:    ", input_path, if (!is.na(input_sheet)) paste0(" [", input_sheet, "]"), "\n", sep = "")
cat("Samples:  ", nrow(results), " included of ", nrow(x), "\n", sep = "")
cat("Age:      ", fmt(min(age), 5), " to ", fmt(max(age), 5), " (", col_age, ")\n", sep = "")
cat("f_GR:     ", if (C33_in_index) "(C31+C33)/(C27+C31+C33)" else "C31/(C27+C31) - no C33 available", "\n", sep = "")
cat("Outputs:  ", normalizePath(output_dir, winslash = "/"), "\n\n", sep = "")
cat("Significant intervals:\n"); cat(interval_lines, sep = "\n")
cat("\nGAM diagnostics:\n")
print(gam_diagnostics[, c("series", "n", "basis_k", "edf", "k_index", "k_check_pvalue",
  "deviance_explained", "residual_lag1_acf", "percent_of_record_significant")], row.names = FALSE, digits = 3)
