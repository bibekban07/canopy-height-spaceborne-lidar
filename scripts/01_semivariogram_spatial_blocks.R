if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, dplyr, sf, gstat, ggplot2)

DATA_DIR    <- file.path(getwd(), "data")
RESULTS_DIR <- file.path(getwd(), "results")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

sample_file <- file.path(DATA_DIR, "sample_training_pts.csv")
if (!file.exists(sample_file)) stop(paste("File not found:", sample_file))

df_pts <- read_csv(sample_file, show_col_types = FALSE)
target_col <- intersect(c("rh98", "target_height_m", "h_canopy"), names(df_pts))[1]
lon_col    <- intersect(c("lon", "longitude", "x"), tolower(names(df_pts)))[1]
lat_col    <- intersect(c("lat", "latitude", "y"), tolower(names(df_pts)))[1]

sf_pts <- df_pts %>%
  filter(!is.na(.data[[target_col]]), !is.na(.data[[lon_col]]), !is.na(.data[[lat_col]])) %>%
  st_as_sf(coords = c(lon_col, lat_col), crs = 4326) %>%
  st_transform(crs = 5070)

form <- as.formula(paste(target_col, "~ 1"))
exp_vgm <- variogram(form, data = as(sf_pts, "Spatial"), cutoff = 10000, width = 500)
init_model <- vgm(psill = var(df_pts[[target_col]], na.rm = TRUE), model = "Sph", range = 4000, nugget = 5)
fitted_vgm <- fit.variogram(exp_vgm, init_model)
print(fitted_vgm)

range_m  <- round(fitted_vgm$range[2], 1)
sill_val <- round(sum(fitted_vgm$psill), 2)
cat(sprintf("
-> Fitted Semivariogram Range: %.1f m (~%.2f km) | Sill: %.2f
", range_m, range_m/1000, sill_val))

plot_path <- file.path(RESULTS_DIR, "Figure_S1_semivariogram_sample.png")
png(plot_path, width = 7, height = 5, units = "in", res = 300)
plot(exp_vgm, model = fitted_vgm, 
     main = sprintf("Empirical Semivariogram (%s)
Range: %.0f m | Sill: %.1f", target_col, range_m, sill_val),
     xlab = "Distance Lag (m)", ylab = "Semivariance (gamma)")
dev.off()
cat(sprintf("-> Plot saved to: %s
", plot_path))
