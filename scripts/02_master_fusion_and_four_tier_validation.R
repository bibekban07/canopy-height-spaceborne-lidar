if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, dplyr, ranger, xgboost)

DATA_DIR    <- file.path(getwd(), "data")
RESULTS_DIR <- file.path(getwd(), "results")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

df_train  <- read_csv(file.path(DATA_DIR, "sample_training_pts.csv"), show_col_types = FALSE)
df_val15  <- if (file.exists(file.path(DATA_DIR, "sample_val15m.csv"))) read_csv(file.path(DATA_DIR, "sample_val15m.csv"), show_col_types = FALSE) else NULL
df_val001 <- if (file.exists(file.path(DATA_DIR, "sample_val001deg.csv"))) read_csv(file.path(DATA_DIR, "sample_val001deg.csv"), show_col_types = FALSE) else NULL

if ("sensor_origin" %in% names(df_train)) {
  df_train$sensor_origin <- as.factor(df_train$sensor_origin)
  if (!is.null(df_val15) && "sensor_origin" %in% names(df_val15)) df_val15$sensor_origin <- factor(df_val15$sensor_origin, levels = levels(df_train$sensor_origin))
  if (!is.null(df_val001) && "sensor_origin" %in% names(df_val001)) df_val001$sensor_origin <- factor(df_val001$sensor_origin, levels = levels(df_train$sensor_origin))
}

get_metrics <- function(actual, predicted, target_name, tier_name, algo) {
  rmse  <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
  r2    <- cor(actual, predicted, use = "complete.obs")^2
  mae   <- mean(abs(actual - predicted), na.rm = TRUE)
  bias  <- mean(predicted - actual, na.rm = TRUE)
  rrmse <- (rmse / mean(actual, na.rm = TRUE)) * 100
  return(data.frame(
    Target = target_name, Algorithm = algo, Validation_Tier = tier_name,
    R2 = round(r2, 4), RMSE = round(rmse, 2), rRMSE = round(rrmse, 2), 
    MAE = round(mae, 2), Bias = round(bias, 2)
  ))
}

preds_top <- c("topo_northness", "topo_eastness", "topo_slope", "topo_elevation",
               "canopy_cover_pct", "sensor_origin", "landsat8_red", "landsat8_green", 
               "landsat8_blue", "landsat8_swir2", "sentinel2_ndre", "sentinel2_rervi", 
               "sentinel2_red", "sentinel1_vh_backscatter", "sentinel1_vv_backscatter")

preds_mean <- c("topo_northness", "topo_eastness", "topo_slope", "topo_elevation",
                "canopy_cover_pct", "sensor_origin", "landsat8_swir2", "landsat8_red", 
                "landsat8_green", "landsat8_blue", "sentinel2_ndre", "sentinel2_red", 
                "sentinel2_swir2", "sentinel1_vh_backscatter", "sentinel1_vv_backscatter")

results_matrix <- data.frame()

for (t_type in c("Max", "Mean")) {
  target_col <- if (t_type == "Max") intersect(c("rh98", "target_height_m", "h_canopy"), names(df_train))[1] else intersect(c("rh50", "target_mean_height_m", "h_mean_canopy"), names(df_train))[1]
  display_target <- if (t_type == "Max") "Top Canopy (H_top)" else "Mean Canopy (H_mean)"
  final_preds    <- if (t_type == "Max") preds_top else preds_mean
  
  # Tier 1: Random Split
  set.seed(42)
  train_idx <- sample(1:nrow(df_train), 0.8 * nrow(df_train))
  d_train   <- df_train[train_idx, ]
  d_test    <- df_train[-train_idx, ]
  
  rf_model  <- ranger(as.formula(paste(target_col, "~", paste(final_preds, collapse = "+"))), data = d_train, num.trees = 500, seed = 42, num.threads = 1)
  xgb_model <- xgb.train(params = list(objective = "reg:squarederror", max_depth = 6, nthread = 1), 
                         data = xgb.DMatrix(data.matrix(d_train[, final_preds]), label = d_train[[target_col]]), nrounds = 100, verbose = 0)
  
  p_rf_t1  <- predict(rf_model, data = d_test)$predictions
  p_xgb_t1 <- predict(xgb_model, xgb.DMatrix(data.matrix(d_test[, final_preds])))
  
  results_matrix <- rbind(results_matrix, 
                          get_metrics(d_test[[target_col]], p_rf_t1, display_target, "Tier_1_Random", "RF"),
                          get_metrics(d_test[[target_col]], p_xgb_t1, display_target, "Tier_1_Random", "XGBoost"))
  
  # Tier 2: Spatial Block (~5 km)
  lon_col <- intersect(c("lon", "longitude", "x"), tolower(names(df_train)))[1]
  lat_col <- intersect(c("lat", "latitude", "y"), tolower(names(df_train)))[1]
  df_spat <- df_train %>% mutate(block_id = as.integer(as.factor(paste(floor(.data[[lon_col]] * 20), floor(.data[[lat_col]] * 20)))))
  
  set.seed(42)
  train_blocks <- sample(unique(df_spat$block_id), 0.8 * length(unique(df_spat$block_id)))
  s_train      <- df_spat %>% filter(block_id %in% train_blocks)
  s_test       <- df_spat %>% filter(!(block_id %in% train_blocks))
  
  rf_spat  <- ranger(as.formula(paste(target_col, "~", paste(final_preds, collapse = "+"))), data = s_train, num.trees = 500, seed = 42, num.threads = 1)
  xgb_spat <- xgb.train(params = list(objective = "reg:squarederror", max_depth = 6, nthread = 1), 
                        data = xgb.DMatrix(data.matrix(s_train[, final_preds]), label = s_train[[target_col]]), nrounds = 100, verbose = 0)
  
  p_rf_t2  <- predict(rf_spat, data = s_test)$predictions
  p_xgb_t2 <- predict(xgb_spat, xgb.DMatrix(data.matrix(s_test[, final_preds])))
  
  results_matrix <- rbind(results_matrix, 
                          get_metrics(s_test[[target_col]], p_rf_t2, display_target, "Tier_2_Spatial_Block", "RF"),
                          get_metrics(s_test[[target_col]], p_xgb_t2, display_target, "Tier_2_Spatial_Block", "XGBoost"))
  
  # Tier 3: 15m Sub-pixel Holdout
  if (!is.null(df_val15)) {
    p_rf_t3  <- predict(rf_model, data = df_val15)$predictions
    p_xgb_t3 <- predict(xgb_model, xgb.DMatrix(data.matrix(df_val15[, final_preds])))
    results_matrix <- rbind(results_matrix, 
                            get_metrics(df_val15[[target_col]], p_rf_t3, display_target, "Tier_3_Holdout_15m", "RF"),
                            get_metrics(df_val15[[target_col]], p_xgb_t3, display_target, "Tier_3_Holdout_15m", "XGBoost"))
  }
  
  # Tier 4: 0.001 deg Stand-scale Holdout
  if (!is.null(df_val001)) {
    p_rf_t4  <- predict(rf_model, data = df_val001)$predictions
    p_xgb_t4 <- predict(xgb_model, xgb.DMatrix(data.matrix(df_val001[, final_preds])))
    results_matrix <- rbind(results_matrix, 
                            get_metrics(df_val001[[target_col]], p_rf_t4, display_target, "Tier_4_Holdout_0.001deg", "RF"),
                            get_metrics(df_val001[[target_col]], p_xgb_t4, display_target, "Tier_4_Holdout_0.001deg", "XGBoost"))
  }
}

print(as.data.frame(results_matrix))
write_csv(results_matrix, file.path(RESULTS_DIR, "Table_S_Validation_Test_Metrics.csv"))
