# =============================================================================
# Echinococcus Covariate Extraction
#
# Purpose:
#   Extract habitat covariates for each Echinococcus sample within the
#   sampled coyote's 95% AKDE home range, or, for animals without a home
#   range, within a circular buffer of equivalent area (see
#   build_echino_buffers.R). All covariates are also z-scaled
#   (mean 0, SD 1; suffix "_scl") for use in models.
#
#   Covariates:
#     prop_<class>  Proportion of each land-cover class (NLCD), from the
#                   annual land-cover layer matching the sample year
#     avg_ndvi      Mean NDVI, from the annual layer matching the sample year
#     avg_ndmi      Mean NDMI, from the annual layer matching the sample year
#     avg_soil      Mean plant-available water capacity, 0-30 cm
#                   (SoilGrids; static)
#
#   All values are area-weighted: each raster cell contributes in proportion
#   to the fraction of it covered by the polygon, so cells on the polygon
#   edge are counted accurately. Land-cover proportions are calculated
#   relative to all non-missing cells in the polygon.
#
# Polygon used for each sample (applied identically to every covariate):
#   1. Home range matching the sample's age class: akde_95_<id><suffix>.shp
#      (suffix: _a = adult, _s = subadult, _p = pup)
#   2. Otherwise, the combined home range:         akde_95_<id>.shp
#   3. Otherwise, the AKDE-size buffer:            buffer_<id><suffix>.shp
#   This matches the rule make_akde_size_buffers.R uses to decide which
#   samples need a buffer, so every sample has exactly one polygon.
#
# Input (paths set in CONFIGURATION below):
#   CSV with one row per sample and at least the columns:
#     id             - animal identifier (matches shapefile names)
#     age_class      - "adult" / "subadult" / "pup" (case-insensitive)
#     date_collected - sample date, containing a 4-digit year
#   Raster layers listed in CONFIGURATION. Annual rasters have one band per
#   year; polygons are reprojected to each raster's CRS before extraction.
#
# Output (written to OUTPUT_CSV; the input file is not modified):
#   All original columns, plus:
#     geom_source - "home_range" or "buffer"
#     geom_file   - name of the shapefile used
#     each covariate above, and its scaled version (<name>_scl)
#
# Requirements:
#   R (>= 4.1), sf, terra
#
# Author:  Laurel Serieys
# Note:    Code was refactored and documented with assistance from
#          Claude (Anthropic), an AI assistant. All analytical decisions,
#          parameter values, and validation are the author's.
# =============================================================================


# ---- 1. SETUP ---------------------------------------------------------------

  library(sf)
  library(terra)


# Paths are relative to the repository root. Edit to match your local layout.
echino_data <- file.path("data", "echinococcus_data.csv")
OUTPUT_CSV  <- file.path("output", "echinococcus_data_covariates.csv")
HR_DIR      <- file.path("data", "akde_home_ranges")
BUFFER_DIR  <- file.path("output", "buffers_akde_size")
RASTER_DIR  <- file.path("data", "rasters")

SUFFIX_LOOKUP <- c(adult = "_a", subadult = "_s", pup = "_p")

# Land cover (annual NLCD; bands named with their year, e.g. "NLCD_2018").
NLCD <- list(
  path  = file.path(RASTER_DIR, "nlcd_annual_2000_2024.tif"),
  years = NULL  # NULL = read each band's year from its name
)

LANDCOVER_CLASSES <- list(
  developed   = c(22, 23, 24),
  open        = c(21),
  forest      = c(41, 42, 43),
  shrub       = c(52),
  agriculture = c(81, 82),
  grassland   = c(71),
  wetland     = c(90, 95),
  water       = c(11)
)

# Annual continuous covariates: the band matching the sample year is used.
# years = NULL reads each band's year from its name; otherwise years are
# assigned to bands in order.
ANNUAL_COVARIATES <- list(
  ndvi = list(path = file.path(RASTER_DIR, "ndvi_annual_2000_2024.tif"), years = NULL),
  ndmi = list(path = file.path(RASTER_DIR, "ndmi_2016_2021.tif"),        years = 2016:2021)
)

# Static continuous covariates: one band (by name or position) is used.
STATIC_COVARIATES <- list(
  soil = list(path = file.path(RASTER_DIR, "soilgrids_vwc_pawc_250m.tif"),   band = "pawc_0_30cm")
)

dir.create(dirname(OUTPUT_CSV), showWarnings = FALSE, recursive = TRUE)


# ---- 2. HELPER FUNCTIONS ----------------------------------------------------

#' Open an annual raster and name its bands by year.
load_annual_raster <- function(path, years = NULL) {
  r <- rast(path)
  if (is.null(years)) {
    years <- as.integer(regmatches(names(r), regexpr("(19|20)[0-9]{2}", names(r))))
    if (length(years) != nlyr(r))
      stop("Could not read a year from every band name in ", path,
           ". Band names: ", paste(names(r), collapse = ", "))
  } else if (length(years) != nlyr(r)) {
    stop(path, " has ", nlyr(r), " bands but ", length(years), " years were given.")
  }
  names(r) <- as.character(years)
  r
}

#' Locate the polygon for one sample: age-matched home range, then combined
#' home range, then AKDE-size buffer. Returns the dissolved geometry.
get_geometry <- function(id, suffix) {
  candidates <- c(
    home_range = file.path(HR_DIR,     paste0("akde_95_", id, suffix, ".shp")),
    home_range = file.path(HR_DIR,     paste0("akde_95_", id, ".shp")),
    buffer     = file.path(BUFFER_DIR, paste0("buffer_",  id, suffix, ".shp"))
  )
  hit <- which(file.exists(candidates))[1]
  if (is.na(hit)) return(NULL)
  
  geom <- st_make_valid(st_union(st_geometry(st_read(candidates[hit], quiet = TRUE))))
  list(geom = geom, source = names(candidates)[hit], file = basename(candidates[hit]))
}

#' Area-weighted cell values within a polygon, for one raster layer.
#' Returns a data frame of cell values and coverage fractions (NA cells removed).
extract_cells <- function(layer, geom) {
  v  <- vect(st_transform(geom, crs(layer)))
  ex <- terra::extract(layer, v, exact = TRUE, ID = FALSE)
  data.frame(value = ex[[1]], fraction = ex$fraction)[!is.na(ex[[1]]), ]
}

#' Area-weighted mean of a raster layer within a polygon.
weighted_mean <- function(layer, geom) {
  ex <- extract_cells(layer, geom)
  if (nrow(ex) == 0 || sum(ex$fraction) == 0) return(NA_real_)
  sum(ex$value * ex$fraction) / sum(ex$fraction)
}

#' Area-weighted proportion of each land-cover class within a polygon.
landcover_props <- function(layer, geom, classes) {
  ex <- extract_cells(layer, geom)
  total <- sum(ex$fraction)
  sapply(classes, function(codes) {
    if (total == 0) return(NA_real_)
    sum(ex$fraction[ex$value %in% codes]) / total
  })
}


# ---- 3. READ AND PREPARE DATA -----------------------------------------------

if (!file.exists(echino_data)) stop("Input file not found: ", echino_data)
if (!dir.exists(HR_DIR))       stop("Home-range directory not found: ", HR_DIR)
if (!dir.exists(BUFFER_DIR))   stop("Buffer directory not found: ", BUFFER_DIR)

df <- read.csv(echino_data, stringsAsFactors = FALSE)

missing_cols <- setdiff(c("id", "age_class", "date_collected"), names(df))
if (length(missing_cols) > 0)
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))

# Shapefile suffix from age class, and sample year from the collection date.
suffix <- unname(SUFFIX_LOOKUP[tolower(trimws(df$age_class))])
if (any(is.na(suffix)))
  stop("Unrecognized age_class for id(s): ", paste(df$id[is.na(suffix)], collapse = ", "))

year_str    <- regmatches(df$date_collected, regexpr("(19|20)[0-9]{2}", df$date_collected))
sample_year <- rep(NA_character_, nrow(df))
sample_year[grepl("(19|20)[0-9]{2}", df$date_collected)] <- year_str
if (any(is.na(sample_year)))
  stop("Could not read a 4-digit year from date_collected for id(s): ",
       paste(df$id[is.na(sample_year)], collapse = ", "))

message(nrow(df), " samples; years ", paste(sort(unique(sample_year)), collapse = ", "))


# ---- 4. LOAD RASTERS --------------------------------------------------------

nlcd_r <- load_annual_raster(NLCD$path, NLCD$years)

annual_r <- lapply(ANNUAL_COVARIATES, function(x) load_annual_raster(x$path, x$years))

static_r <- lapply(STATIC_COVARIATES, function(x) {
  r <- rast(x$path)
  if (is.character(x$band) && !x$band %in% names(r))
    stop("Band '", x$band, "' not found in ", x$path,
         ". Available bands: ", paste(names(r), collapse = ", "))
  r[[x$band]]
})

# Every sample year must be available in every annual raster.
for (nm in c("nlcd", names(annual_r))) {
  r <- if (nm == "nlcd") nlcd_r else annual_r[[nm]]
  absent <- setdiff(unique(sample_year), names(r))
  if (length(absent) > 0)
    stop("Sample year(s) ", paste(absent, collapse = ", "), " not in the ", nm, " raster.")
}


# ---- 5. EXTRACT COVARIATES --------------------------------------------------

lc_cols     <- paste0("prop_", names(LANDCOVER_CLASSES))
cont_cols   <- paste0("avg_", c(names(ANNUAL_COVARIATES), names(STATIC_COVARIATES)))
covariates  <- matrix(NA_real_, nrow = nrow(df), ncol = length(lc_cols) + length(cont_cols),
                      dimnames = list(NULL, c(lc_cols, cont_cols)))
geom_source <- rep(NA_character_, nrow(df))
geom_file   <- rep(NA_character_, nrow(df))

for (i in seq_len(nrow(df))) {
  g <- get_geometry(df$id[i], suffix[i])
  
  if (is.null(g)) {
    warning("No home range or buffer found for id ", df$id[i], suffix[i], "; covariates set to NA.")
    next
  }
  geom_source[i] <- g$source
  geom_file[i]   <- g$file
  yr <- sample_year[i]
  
  covariates[i, lc_cols] <- landcover_props(nlcd_r[[yr]], g$geom, LANDCOVER_CLASSES)
  
  for (nm in names(annual_r)) {
    covariates[i, paste0("avg_", nm)] <- weighted_mean(annual_r[[nm]][[yr]], g$geom)
  }
  for (nm in names(static_r)) {
    covariates[i, paste0("avg_", nm)] <- weighted_mean(static_r[[nm]], g$geom)
  }
  
  # All-NA values usually mean the polygon does not overlap the rasters
  # (e.g. a shapefile with an incorrect CRS).
  if (all(is.na(covariates[i, ])))
    warning("No raster values found within ", g$file, "; check its CRS.")
  
  message(sprintf("[%d/%d] id %s%s  %s  (%s)", i, nrow(df), df$id[i], suffix[i],
                  g$file, yr))
}


# ---- 6. SCALE AND SAVE ------------------------------------------------------
# Covariates are z-scaled across samples; missing values stay NA.

df_out <- cbind(df, geom_source = geom_source, geom_file = geom_file,
                as.data.frame(covariates))

for (col in colnames(covariates)) {
  df_out[[paste0(col, "_scl")]] <- as.numeric(scale(df_out[[col]]))
}

write.csv(df_out, OUTPUT_CSV, row.names = FALSE)

message("\nDone. Covariates for ", sum(!is.na(geom_source)), " of ", nrow(df),
        " samples written to ", OUTPUT_CSV)
print(table(geom_source, useNA = "ifany"))
print(summary(df_out[, colnames(covariates)]))


# ---- 7. SESSION INFO --------------------------------------------------------
# Record package versions for reproducibility.

sessionInfo()