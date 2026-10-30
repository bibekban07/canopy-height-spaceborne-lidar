if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, dplyr, ranger, ggplot2)

DATA_DIR    <- file.path(getwd(), "data")
RESULTS_DIR <- file.path(getwd(), "results")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

train_file  <- file.path(DATA_DIR, "sample_training_pts.csv")
als_pt_file <- file.path(DATA_DIR, "sample_als_validation_pts.csv")

stopifnot(file.exists(train_file), file.exists(als_pt_file))

df_train <- read_csv(train_file, show_col_types = FALSE)
df_als   <- read_csv(als_pt_file, show_col_types = FALSE)

if ("sensor_origin" %in% names(df_train)) df_train$sensor_origin <- as.factor(df_train$sensor_origin)
if ("sensor_origin" %in% names(df_als)) df_als$sensor_origin <- factor(df_als$sensor_origin, levels = levels(df_train$sensor_origin))

preds_top <- c("topo_northness", "topo_eastness", "topo_slope", "topo_elevation",
               "canopy_cover_pct", "sensor_origin", "landsat8_red", "landsat8_green", 
               "landsat8_blue", "landsat8_swir2", "sentinel2_ndre", "sentinel2_rervi", 
               "sentinel2_red", "sentinel1_vh_backscatter", "sentinel1_vv_backscatter")

preds_mean <- c("topo_northness", "topo_eastness", "topo_slope", "topo_elevation",
                "canopy_cover_pct", "sensor_origin", "landsat8_swir2", "landsat8_red", 
                "landsat8_green", "landsat8_blue", "sentinel2_ndre", "sentinel2_red", 
                "sentinel2_swir2", "sentinel1_vh_backscatter", "sentinel1_vv_backscatter")

get_metrics <- function(actual, predicted, eval_type, target_name) {
  rmse  <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
  r2    <- cor(actual, predicted, use = "complete.obs")^2
  mae   <- mean(abs(actual - predicted), na.rm = TRUE)
  bias  <- mean(predicted - actual, na.rm = TRUE)
  rrmse <- (rmse / mean(actual, na.rm = TRUE)) * 100
  return(data.frame(
    Evaluation_Type = eval_type, Target = target_name, N_Points = length(na.omit(actual)),
    R2 = round(r2, 4), RMSE = round(rmse, 2), rRMSE = round(rrmse, 2), MAE = round(mae, 2), Bias = round(bias, 2),
    Remarks = "ALS Point Track Holdout (Wall-to-wall raster grids omitted due to repository storage limits)"
  ))
}

compiled_als <- data.frame()

for (tgt in c("Top", "Mean")) {
  target_col   <- if (tgt == "Top") intersect(c("rh98", "target_height_m", "h_canopy"), names(df_train))[1] else intersect(c("rh50", "target_mean_height_m", "h_mean_canopy"), names(df_train))[1]
  final_preds  <- if (tgt == "Top") preds_top else preds_mean
  clean_target <- if (tgt == "Top") "Top Canopy (H_top)" else "Mean Canopy (H_mean)"
  obs_col      <- if (tgt == "Top") intersect(c("Observed_ALS_Top", "Observed_ALS"), names(df_als))[1] else intersect(c("Observed_ALS_Mean", "Observed_ALS"), names(df_als))[1]
  
  set.seed(42)
  rf_model <- ranger(as.formula(paste(target_col, "~", paste(final_preds, collapse = "+"))),
                     data = df_train, num.trees = 500, seed = 42, num.threads = 1)
  
  pred_als  <- predict(rf_model, data = df_als)$predictions
  obs_als   <- df_als[[obs_col]]
  valid_idx <- which(!is.na(obs_als) & !is.na(pred_als) & obs_als > 0)
  
  m <- get_metrics(obs_als[valid_idx], pred_als[valid_idx], "ALS_Point_Holdout", clean_target)
  compiled_als <- rbind(compiled_als, m)
  
  max_val <- if (tgt == "Top") 48 else 35
  breaks_seq <- if (tgt == "Top") seq(0, 40, by = 10) else seq(0, 30, by = 10)
  
  p <- ggplot(data.frame(Observed = obs_als[valid_idx], Estimated = pred_als[valid_idx]), aes(x = Observed, y = Estimated)) +
    geom_point(shape = 1, alpha = 0.35, color = "black", size = 1.6) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 0.9) +
    annotate("text", x = 2, y = max_val - 6, hjust = 0, size = 4.2,
             label = sprintf("R² = %.2f
RMSE = %.2f m
rRMSE = %.2f%%
MAE = %.2f m
Bias = %.2f m
(N = %d points)", 
                             m$R2, m$RMSE, m$rRMSE, m$MAE, m$Bias, m$N_Points)) +
    scale_x_continuous(limits = c(0, max_val), breaks = breaks_seq) +
    scale_y_continuous(limits = c(0, max_val), breaks = breaks_seq) +
    theme_bw() + 
    labs(title = sprintf("Independent ALS Point Validation (30m %s)", clean_target),
         subtitle = "Trained on 5,000 Demonstration Sample | Evaluated on ALS Overlap Tracks",
         x = sprintf("Observed ALS Reference Height (m) [%s]", obs_col), y = "Estimated Height (m)") + 
    coord_fixed(ratio = 1)
  
  ggsave(file.path(RESULTS_DIR, sprintf("Sample_ALS_Point_Validation_%s.png", tgt)), plot = p, width = 6.5, height = 6.5, dpi = 300)
}

print(as.data.frame(compiled_als))
write_csv(compiled_als, file.path(RESULTS_DIR, "Sample_ALS_Point_Validation_Table.csv"))
