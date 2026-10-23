# =============================================================================
# 95% AKDE home range estimation for Chicago coyotes (amt + ctmm)
#
# For each individual, this script:
#   1. Fits continuous-time movement models with ctmm and selects the best
#      model by AICc (ctmm.select).
#   2. Estimates the 95% autocorrelated kernel density estimate (AKDE) home
#      range from the selected model.
#   3. Clips the 95% isopleth to remove area over Lake Michigan.
#   4. Writes one shapefile per animal, a summary table, and an .rds file of
#      all fitted AKDE objects.
#
# Coordinates: input in UTM Zone 16N NAD27 (EPSG:26716), reprojected to
#              UTM Zone 16N WGS84 (EPSG:32616) for all analysis and outputs.
# Timestamps:  input is local Chicago time; parsed directly in
#              America/Chicago. Times that do not exist (DST spring-forward
#              gap) fail to parse and are dropped.
# IDs:         optionally pool pup (_p) fixes with the same animal's subadult
#              (_s) fixes (see `merge_pup_subadult` below; off by default and
#              off for the published analysis).
#
# Author:  Laurel Serieys
# Note:    Code was refactored and documented with assistance from
#          Claude (Anthropic), an AI assistant. All analytical decisions,
#          parameter values, and validation are the author's.
# License: see LICENSE file in this repository
# =============================================================================

  library(amt)
  library(ctmm)
  library(sf)
  library(terra)
  library(dplyr)
  library(lubridate)
  library(readr)

# -----------------------------------------------------------------------------
# CONFIGURATION — edit these to match your setup
# -----------------------------------------------------------------------------

# Input files (paths relative to the repository root)
data_path      <- "data/GPS_VHF_Locations.csv"     # columns: id, easting, northing, datetime
lake_mask_path <- "data/lake_michigan_mask.tif"     # raster; 1 = lake, other = land

# Output directory (created if it does not exist)
output_dir     <- "output"

# Coordinate reference systems
input_crs <- 26716  # Input data: UTM Zone 16N NAD27
utm_crs   <- 32616  # Working/output CRS: UTM Zone 16N WGS84

# Home range settings
ud_level      <- 0.95  # isopleth level of the utilization distribution
min_fixes     <- 20    # animals with fewer fixes are skipped
dof_threshold <- 15    # animals with effective sample size (DOF[area]) below
                       # this are still estimated but flagged

# Pool _p (pup) fixes into the matching _s (subadult) ID before estimation?
merge_pup_subadult <- FALSE

# -----------------------------------------------------------------------------
# SETUP
# -----------------------------------------------------------------------------

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}

summary_csv <- file.path(output_dir, "AKDE_homerange_summary.csv")
rds_file    <- file.path(output_dir, "all_akdes.rds")
session_txt <- file.path(output_dir, "sessionInfo.txt")

# -----------------------------------------------------------------------------
# LOAD & PREPARE LAKE MICHIGAN MASK
# -----------------------------------------------------------------------------

message("Loading Lake Michigan mask from: ", lake_mask_path)
lake_rast      <- rast(lake_mask_path)
lake_rast_proj <- project(lake_rast, paste0("EPSG:", utm_crs))

lake_poly <- as.polygons(lake_rast_proj == 1, dissolve = TRUE) %>%
  st_as_sf() %>%
  filter(.[[1]] == 1) %>%
  st_union() %>%
  st_make_valid()

message("  Lake polygon CRS: EPSG:", st_crs(lake_poly)$epsg)

# -----------------------------------------------------------------------------
# LOAD & PREPARE LOCATION DATA
# -----------------------------------------------------------------------------

message("\nLoading location data from: ", data_path)
raw <- read.csv(data_path, stringsAsFactors = FALSE)

message("  Columns found: ", paste(names(raw), collapse = ", "))
message("  Total rows: ", nrow(raw))

required_cols <- c("id", "easting", "northing", "datetime")
missing_cols  <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}

raw <- raw %>% rename(x = easting, y = northing)

# Timestamps are local Chicago time. Times falling in the DST spring-forward
# gap do not exist in America/Chicago, become NA, and are dropped below.
raw$datetime <- parse_date_time(
  raw$datetime,
  orders = c("ymd HMS", "ymd HM", "mdy HMS", "mdy HM", "ymd"),
  tz     = "America/Chicago"
)

if (all(is.na(raw$datetime))) {
  stop("Column 'datetime' could not be parsed. Check the timestamp format in the input CSV.")
}

n_failed <- sum(is.na(raw$datetime))
if (n_failed > 0) {
  message("  WARNING: ", n_failed,
          " timestamp(s) failed to parse (e.g. DST gap) and will be dropped.")
}

# Shift fixes recorded at exactly 00:00:00 forward by one minute to prevent parsing issues.
midnight   <- format(raw$datetime, "%H:%M:%S") == "00:00:00"
n_midnight <- sum(midnight, na.rm = TRUE)
if (n_midnight > 0) {
  raw$datetime[which(midnight)] <- raw$datetime[which(midnight)] + 60
  message("  Adjusted ", n_midnight, " midnight timestamp(s) from 00:00 to 00:01.")
}

# Drop rows with missing location or time
n_before <- nrow(raw)
raw <- raw %>% filter(!is.na(x), !is.na(y), !is.na(datetime))
n_dropped <- n_before - nrow(raw)
if (n_dropped > 0) {
  message("  Dropped ", n_dropped, " row(s) with missing coordinates or timestamps.")
}

# -----------------------------------------------------------------------------
# OPTIONAL: COMBINE _p (PUP) AND _s (SUBADULT) RECORDS FOR THE SAME ANIMAL
# For any animal with both <id>_p and <id>_s, the _p rows are relabeled _s.
# -----------------------------------------------------------------------------

if (merge_pup_subadult) {
  all_ids       <- unique(raw$id)
  pup_base      <- sub("_p$", "", all_ids[grepl("_p$", all_ids)])
  sub_base      <- sub("_s$", "", all_ids[grepl("_s$", all_ids)])
  matching_base <- intersect(pup_base, sub_base)

  if (length(matching_base) == 0) {
    message("\nNo animals with both _p and _s records; no IDs merged.")
  } else {
    message("\nMerging _p and _s records for ", length(matching_base), " animal(s):")
    for (base in matching_base) {
      id_p <- paste0(base, "_p")
      id_s <- paste0(base, "_s")
      n_p  <- sum(raw$id == id_p)
      n_s  <- sum(raw$id == id_s)
      raw$id[raw$id == id_p] <- id_s
      message("  ", id_p, " (", n_p, ") + ", id_s, " (", n_s, ") -> ",
              id_s, " (", n_p + n_s, " fixes)")
    }
  }
}

# -----------------------------------------------------------------------------
# REPROJECT EPSG:26716 -> EPSG:32616
# -----------------------------------------------------------------------------

message("\nReprojecting coordinates from EPSG:", input_crs, " to EPSG:", utm_crs, "...")

pts_proj <- st_as_sf(raw, coords = c("x", "y"), crs = input_crs) %>%
  st_transform(crs = utm_crs)

xy    <- st_coordinates(pts_proj)
raw$x <- xy[, 1]
raw$y <- xy[, 2]

animal_ids <- sort(unique(raw$id))
message("\nFound ", length(animal_ids), " individual(s): ",
        paste(animal_ids, collapse = ", "))

# -----------------------------------------------------------------------------
# HELPER: empty summary row for skipped / failed animals
# -----------------------------------------------------------------------------

empty_summary <- function(aid, status, n_fixes = NA_integer_, duration_days = NA_real_) {
  data.frame(
    id               = as.character(aid),
    model_selected   = NA_character_,
    n_fixes          = n_fixes,
    duration_days    = duration_days,
    tau_pos_months   = NA_real_,
    dof_area         = NA_real_,
    area_low_km2     = NA_real_,
    area_est_km2     = NA_real_,
    area_high_km2    = NA_real_,
    area_clipped_km2 = NA_real_,
    dof_flag         = NA_character_,
    shapefile        = NA_character_,
    status           = status,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# LOOP — AKDE HOME RANGE PER INDIVIDUAL
# -----------------------------------------------------------------------------

results_list <- lapply(seq_along(animal_ids), function(i) {

  aid <- animal_ids[i]
  message("\n--- Processing: ", aid, " (", i, " of ", length(animal_ids), ") ---")

  tryCatch({

    ind_data <- raw %>%
      filter(id == aid) %>%
      arrange(datetime)

    n_fixes       <- nrow(ind_data)
    duration_days <- round(as.numeric(difftime(max(ind_data$datetime),
                                               min(ind_data$datetime),
                                               units = "days")), 1)

    message("  Fixes: ", n_fixes, " | Duration: ", duration_days, " days")

    if (n_fixes < min_fixes) {
      message("  SKIPPED: fewer than ", min_fixes, " fixes.")
      return(list(
        summary = empty_summary(aid, "SKIPPED - too few fixes", n_fixes, duration_days),
        akde    = NULL
      ))
    }

    # amt track in EPSG:32616
    trk <- make_track(ind_data, .x = x, .y = y, .t = datetime, crs = utm_crs)

    # ── Fit movement models and select best by AICc ─────────────────────────
    message("  Fitting ctmm models...")
    tel        <- as_telemetry(trk)
    guess      <- ctmm.guess(tel, interactive = FALSE)
    m          <- ctmm.select(tel, guess, verbose = FALSE)
    model_name <- summary(m)$name

    # tau[position] in months (seconds / seconds-per-month); NA for IID models
    tau_pos_months <- tryCatch({
      tau_val <- m$tau
      if (length(tau_val) == 0) NA_real_
      else as.numeric(tau_val[1]) / (60 * 60 * 24 * 30.44)
    }, error = function(e) NA_real_)

    # ── AKDE ────────────────────────────────────────────────────────────────
    message("  Computing ", ud_level * 100, "% AKDE...")
    akde_hr <- hr_akde(trk, model = m, level_ud = ud_level)

    hr_poly <- hr_isopleths(akde_hr) %>%
      filter(what == "estimate") %>%
      st_set_crs(utm_crs) %>%
      st_make_valid()

    # Area point estimate and 95% CI from ctmm (pre-clip). units = FALSE
    # returns SI units (m^2) so the conversion to km^2 is unambiguous.
    ctmm_summary  <- summary(akde_hr$akde, level.UD = ud_level, units = FALSE)
    dof_area      <- round(as.numeric(ctmm_summary$DOF["area"]), 2)
    area_ci       <- ctmm_summary$CI[grep("area", rownames(ctmm_summary$CI),
                                          ignore.case = TRUE), , drop = FALSE]
    area_low_km2  <- round(area_ci[1, "low"]  / 1e6, 4)
    area_est_km2  <- round(area_ci[1, "est"]  / 1e6, 4)
    area_high_km2 <- round(area_ci[1, "high"] / 1e6, 4)

    dof_flag <- if (is.na(dof_area) || dof_area < dof_threshold) {
      paste0("LOW DOF (<", dof_threshold, ")")
    } else {
      "OK"
    }

    message("  Model: ", model_name,
            " | DOF(area): ", dof_area, " [", dof_flag, "]",
            " | tau[position]: ",
            ifelse(is.na(tau_pos_months), "N/A",
                   paste0(round(tau_pos_months, 2), " months")))
    message("  AKDE area (pre-clip): ", area_est_km2, " km2 [",
            area_low_km2, " - ", area_high_km2, "]")

    # ── Clip to exclude Lake Michigan ───────────────────────────────────────
    hr_clipped <- tryCatch({
      st_difference(hr_poly, lake_poly) %>% st_make_valid()
    }, error = function(e) {
      message("  WARNING: lake clipping failed (", conditionMessage(e),
              "); using unclipped polygon.")
      hr_poly
    })

    area_clipped_km2 <- round(sum(as.numeric(st_area(hr_clipped))) / 1e6, 4)
    message("  AKDE area (post-clip): ", area_clipped_km2, " km2")

    # ── Attach metadata and write shapefile ─────────────────────────────────
    hr_clipped <- hr_clipped %>%
      mutate(
        id               = as.character(aid),
        level            = ud_level,
        n_fixes          = n_fixes,
        model_selected   = model_name,
        dof_area         = dof_area,
        area_est_km2     = area_est_km2,
        area_clipped_km2 = area_clipped_km2,
        dof_flag         = dof_flag
      ) %>%
      select(id, level, n_fixes, model_selected, dof_area,
             area_est_km2, area_clipped_km2, dof_flag, geometry)

    safe_id <- gsub("[^A-Za-z0-9_-]", "_", as.character(aid))
    out_shp <- file.path(output_dir, paste0("akde_", ud_level * 100, "_", safe_id, ".shp"))
    st_write(hr_clipped, dsn = out_shp, delete_dsn = TRUE, quiet = TRUE)
    message("  Saved: ", basename(out_shp))

    list(
      summary = data.frame(
        id               = as.character(aid),
        model_selected   = model_name,
        n_fixes          = n_fixes,
        duration_days    = duration_days,
        tau_pos_months   = round(tau_pos_months, 2),
        dof_area         = dof_area,
        area_low_km2     = area_low_km2,
        area_est_km2     = area_est_km2,
        area_high_km2    = area_high_km2,
        area_clipped_km2 = area_clipped_km2,
        dof_flag         = dof_flag,
        shapefile        = basename(out_shp),
        status           = "OK",
        stringsAsFactors = FALSE
      ),
      akde = akde_hr
    )

  }, error = function(e) {
    message("  ERROR: ", conditionMessage(e))
    list(summary = empty_summary(aid, paste0("ERROR - ", conditionMessage(e))),
         akde    = NULL)
  })
})

# -----------------------------------------------------------------------------
# SAVE OUTPUTS
# -----------------------------------------------------------------------------

# 1. Summary table (one row per animal)
results_df <- bind_rows(lapply(results_list, `[[`, "summary")) %>%
  arrange(id) %>%
  as_tibble()

write_csv(results_df, summary_csv)
message("\nSummary table saved: ", summary_csv)

# 2. All fitted AKDE objects in a single named list
akde_objects <- setNames(lapply(results_list, `[[`, "akde"),
                         vapply(results_list, function(r) r$summary$id, character(1)))
akde_objects <- Filter(Negate(is.null), akde_objects)

saveRDS(akde_objects, rds_file)
message("AKDE objects saved: ", rds_file, " (", length(akde_objects), " animals)")

# 3. R session info for reproducibility
writeLines(capture.output(sessionInfo()), session_txt)
message("Session info saved: ", session_txt)

# -----------------------------------------------------------------------------
# CONSOLE SUMMARY
# -----------------------------------------------------------------------------

message("\n========== PROCESSING COMPLETE ==========")
print(results_df, n = Inf)

low_dof <- results_df %>% filter(dof_flag == paste0("LOW DOF (<", dof_threshold, ")"))
if (nrow(low_dof) > 0) {
  message("\nAnimals with DOF < ", dof_threshold, " (still included in all outputs):")
  print(low_dof %>% select(id, model_selected, n_fixes, dof_area,
                           area_est_km2, area_clipped_km2), n = Inf)
} else {
  message("\nAll estimated animals have DOF >= ", dof_threshold, ".")
}
