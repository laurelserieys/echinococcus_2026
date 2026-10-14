# =============================================================================
# make_akde_size_buffers.R
#
# Purpose:
#   For coyote Echinococcus sampling records that lack an individual AKDE
#   (autocorrelated kernel density estimate) home range, generate a circular
#   buffer around the sample location whose land area equals the mean 95% AKDE
#   home-range area for that animal's age class and sex. Where a buffer
#   overlaps Lake Michigan, the radius is expanded so that the land-only area
#   (after clipping out the lake) still matches the target area.
#
# Inputs (paths set in the CONFIGURATION section below):
#   1. Echinococcus data CSV with, at minimum, the columns:
#        id         - unique animal identifier
#        sex        - "male" / "female" (case-insensitive)
#        age_class  - "adult" / "subadult" / "pup" (case-insensitive)
#        easting    - x coordinate in the projected CRS given by CRS_EPSG
#        northing   - y coordinate in the projected CRS given by CRS_EPSG
#   2. Directory of existing AKDE home-range shapefiles, named
#        akde_95_<id><suffix>.shp  or  akde_95_<id>.shp
#      where <suffix> is "_a" (adult), "_s" (subadult), or "_p" (pup).
#   3. Lake mask raster (GeoTIFF) in the same CRS, coded 1 = lake, 0 = land.
#
# Outputs (written to OUTPUT_DIR):
#   - buffer_<id><suffix>.shp : one buffer polygon per record without an AKDE
#   - buffer_summary.csv      : target area, achieved area, and radius per buffer
#
# Requirements:
#   R (>= 4.1), sf, terra, dplyr
#
# Author:  Laurel Serieys
# Note:    Code was refactored and documented with assistance from
#          Claude (Anthropic), an AI assistant. All analytical decisions,
#          parameter values, and validation are the author's.
# =============================================================================

  library(sf)
  library(terra)
  library(dplyr)

# ---- CONFIGURATION ----------------------------------------------------------
# Paths are relative to the repository root. Edit to match your local layout.

echino_data <- file.path("data", "echinococcus_data.csv")
HR_DIR      <- file.path("data", "akde_home_ranges")
LAKE_TIF    <- file.path("data", "lake_michigan_mask.tif")
OUTPUT_DIR  <- file.path("output", "buffers_akde_size")

CRS_EPSG <- 26716  # NAD27 / UTM zone 16N

# Target buffer areas (km^2): mean 95% AKDE home-range size by age class and sex.
BUFFER_LOOKUP <- data.frame(
  age_class = c("adult", "adult",  "subadult", "subadult", "pup",  "pup"),
  sex       = c("male",  "female", "male",     "female",   "male", "female"),
  area_km2  = c(11.18,   10.33,    11.10,      16.59,      8.08,   11.49),
  stringsAsFactors = FALSE
)

# File-name suffix used for each age class.
SUFFIX_LOOKUP <- c(adult = "_a", subadult = "_s", pup = "_p")

# Root-finding settings.
MIN_SEARCH_WINDOW_M  <- 3000  # minimum half-width (m) of the lake search window
MAX_WINDOW_DOUBLINGS <- 8     # max times the search window may be doubled
RADIUS_TOL_M         <- 1     # tolerance (m) for the buffer radius solution

# ---- HELPER FUNCTIONS -------------------------------------------------------

#' Return lake polygons (sf) within a square window centered on (x, y),
#' or NULL if no lake cells fall inside the window.
lake_poly_near <- function(lake_r, x, y, half_win) {
  win <- ext(x - half_win, x + half_win, y - half_win, y + half_win)
  win <- terra::intersect(win, ext(lake_r))
  if (is.null(win)) return(NULL)
  
  r_crop <- crop(lake_r, win)
  r_crop[r_crop == 0] <- NA
  if (all(is.na(values(r_crop)))) return(NULL)
  
  st_as_sf(as.polygons(r_crop, dissolve = TRUE))
}

#' Remove lake area from a buffer polygon.
clip_to_land <- function(buf, lake_sf) {
  if (!is.null(lake_sf) && nrow(lake_sf) > 0) {
    overlaps <- suppressWarnings(st_intersects(buf, lake_sf, sparse = FALSE))
    if (any(overlaps)) {
      buf <- suppressWarnings(st_difference(buf, st_union(st_geometry(lake_sf))))
    }
  }
  st_make_valid(buf)
}

#' Land-only area (km^2) of a circular buffer of the given radius (m).
land_area_km2 <- function(pt, radius, lake_sf) {
  buf <- clip_to_land(st_buffer(pt, dist = radius), lake_sf)
  as.numeric(sum(st_area(buf))) / 1e6
}

#' Find the buffer radius (m) whose land-only area equals target_km2.
#' Returns a list with the radius and the lake polygons used.
solve_radius <- function(pt, x, y, target_km2, lake_r) {
  r0  <- sqrt(target_km2 * 1e6 / pi)  # radius if there were no lake overlap
  win <- max(MIN_SEARCH_WINDOW_M, r0 * 3)
  doublings <- 0
  
  repeat {
    lake_sf <- lake_poly_near(lake_r, x, y, win)
    
    # No lake nearby, or no lake overlap at r0: the plain circle is correct.
    if (is.null(lake_sf) || land_area_km2(pt, r0, lake_sf) >= target_km2) {
      return(list(radius = r0, lake_sf = lake_sf))
    }
    
    f <- function(r) land_area_km2(pt, r, lake_sf) - target_km2
    
    # If even a radius equal to the window can't reach the target, widen it.
    if (f(win) < 0 && doublings < MAX_WINDOW_DOUBLINGS) {
      win <- win * 2
      doublings <- doublings + 1
      next
    }
    
    radius <- tryCatch(
      uniroot(f, lower = r0, upper = win, tol = RADIUS_TOL_M)$root,
      error = function(e) {
        warning(sprintf("Root-finding failed; using window radius %.0f m.", win))
        win
      }
    )
    return(list(radius = radius, lake_sf = lake_sf))
  }
}

# ---- INPUT CHECKS -----------------------------------------------------------

for (p in c(echino_data, LAKE_TIF)) {
  if (!file.exists(p)) stop("Input file not found: ", p)
}
if (!dir.exists(HR_DIR)) stop("Home-range directory not found: ", HR_DIR)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- READ AND PREPARE RECORDS ----------------------------------------------

echino_df <- read.csv(echino_data, stringsAsFactors = FALSE)

required_cols <- c("id", "sex", "age_class", "easting", "northing")
missing_cols  <- setdiff(required_cols, names(echino_df))
if (length(missing_cols) > 0) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}

echino_df <- echino_df %>%
  mutate(
    age_class = tolower(trimws(age_class)),
    sex       = tolower(trimws(sex))
  ) %>%
  left_join(BUFFER_LOOKUP, by = c("age_class", "sex")) %>%
  mutate(suffix = unname(SUFFIX_LOOKUP[age_class]))

bad <- is.na(echino_df$area_km2) | is.na(echino_df$suffix)
if (any(bad)) {
  stop("Unrecognized age_class/sex combination for id(s): ",
       paste(echino_df$id[bad], collapse = ", "))
}

# Flag records that already have an AKDE home range (with or without suffix).
echino_df <- echino_df %>%
  mutate(
    has_hr = file.exists(file.path(HR_DIR, paste0("akde_95_", id, suffix, ".shp"))) |
      file.exists(file.path(HR_DIR, paste0("akde_95_", id, ".shp"))),
    needs_buffer = !has_hr
  )

to_buffer <- filter(echino_df, needs_buffer)
message(nrow(to_buffer), " of ", nrow(echino_df),
        " records need a buffer (no matching AKDE home range found).")

# ---- BUILD BUFFERS ----------------------------------------------------------

lake_r  <- rast(LAKE_TIF)
results <- vector("list", nrow(to_buffer))

for (i in seq_len(nrow(to_buffer))) {
  rec <- to_buffer[i, ]
  pt  <- st_sfc(st_point(c(rec$easting, rec$northing)), crs = CRS_EPSG)
  r0  <- sqrt(rec$area_km2 * 1e6 / pi)
  
  sol <- solve_radius(pt, rec$easting, rec$northing, rec$area_km2, lake_r)
  
  buf_final    <- clip_to_land(st_buffer(pt, dist = sol$radius), sol$lake_sf)
  final_km2    <- as.numeric(sum(st_area(buf_final))) / 1e6
  lake_clipped <- !is.null(sol$lake_sf) && sol$radius > r0
  
  out_sf <- st_sf(
    id         = rec$id,
    sex        = rec$sex,
    age_class  = rec$age_class,
    target_km2 = rec$area_km2,
    final_km2  = final_km2,
    radius_m   = sol$radius,
    geometry   = st_geometry(buf_final),
    crs        = CRS_EPSG
  )
  
  out_file <- file.path(OUTPUT_DIR, paste0("buffer_", rec$id, rec$suffix, ".shp"))
  suppressWarnings(st_write(out_sf, out_file, quiet = TRUE, append = FALSE))
  
  results[[i]] <- data.frame(
    id               = rec$id,
    suffix           = rec$suffix,
    sex              = rec$sex,
    age_class        = rec$age_class,
    target_km2       = rec$area_km2,
    final_km2        = final_km2,
    radius_m         = sol$radius,
    clipped_for_lake = lake_clipped
  )
  
  message(sprintf("[%d/%d] id %s%s  target = %.2f km2  final = %.2f km2  r = %.0f m%s",
                  i, nrow(to_buffer), rec$id, rec$suffix, rec$area_km2,
                  final_km2, sol$radius,
                  if (lake_clipped) "  (lake-clipped)" else ""))
}

# ---- SUMMARY ----------------------------------------------------------------

summary_df <- bind_rows(results)
write.csv(summary_df, file.path(OUTPUT_DIR, "buffer_summary.csv"), row.names = FALSE)

message("\nDone. Wrote ", nrow(summary_df), " buffer shapefiles and buffer_summary.csv to ",
        OUTPUT_DIR)
message("Records expanded for lake overlap: ", sum(summary_df$clipped_for_lake))