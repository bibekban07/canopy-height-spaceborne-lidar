if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, dplyr, ranger, ggplot2)

DATA_DIR    <- file.path(getwd(), "data")
RESULTS_DIR <- file.path(getwd(), "results")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

df_train <- read_csv(file.path(DATA_DIR, "sample_training_pts.csv"), show_col_types = FALSE)
if ("sensor_origin" %in% names(df_train)) df_train$sensor_origin <- as.factor(df_train$sensor_origin)

preds_top <- c("topo_northness", "topo_eastness", "topo_slope", "topo_elevation",
               "canopy_cover_pct", "sensor_origin", "landsat8_red", "landsat8_green", 
               "landsat8_blue", "landsat8_swir2", "sentinel2_ndre", "sentinel2_rervi", 
               "sentinel2_red", "sentinel1_vh_backscatter", "sentinel1_vv_backscatter")

preds_mean <- c("topo_northness", "topo_eastness", "topo_slope", "topo_elevation",
                "canopy_cover_pct", "sensor_origin", "landsat8_swir2", "landsat8_red", 
                "landsat8_green", "landsat8_blue", "sentinel2_ndre", "sentinel2_red", 
                "sentinel2_swir2", "sentinel1_vh_backscatter", "sentinel1_vv_backscatter")

for (tgt in c("Top", "Mean")) {
  target_col <- if (tgt == "Top") intersect(c("rh98", "target_height_m", "h_canopy"), names(df_train))[1] else intersect(c("rh50", "target_mean_height_m", "h_mean_canopy"), names(df_train))[1]
  display_title <- if (tgt == "Top") "30m Top Canopy" else "30m Mean Canopy"
  current_preds <- if (tgt == "Top") preds_top else preds_mean
  
  set.seed(42)
  rf_pfi <- ranger(
    formula     = as.formula(paste(target_col, "~", paste(current_preds, collapse = "+"))),
    data        = df_train, num.trees = 500, importance = "permutation", seed = 42, num.threads = 1
  )
  
  imp_df <- data.frame(Variable = names(rf_pfi$variable.importance), Importance = round(as.numeric(rf_pfi$variable.importance), 4))
  imp_df$Category <- "Baseline Control"
  imp_df$Category[grepl("landsat", imp_df$Variable, ignore.case = TRUE)]      <- "Optical (Landsat)"
  imp_df$Category[grepl("sentinel2", imp_df$Variable, ignore.case = TRUE)]    <- "Red-Edge (Sentinel-2)"
  imp_df$Category[grepl("sentinel1ert{}sarert{}vvert{}vh", imp_df$Variable, ignore.case = TRUE)] <- "Radar (Sentinel-1)"
  imp_df$Category[grepl("topo", imp_df$Variable, ignore.case = TRUE)]         <- "Topography"
  
  imp_plot <- imp_df %>% arrange(Importance)
  imp_plot$Variable <- factor(imp_plot$Variable, levels = imp_plot$Variable)
  
  p_color <- ggplot(imp_plot, aes(x = Variable, y = Importance, fill = Category)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    coord_flip() + theme_bw() +
    labs(title = sprintf("Permutation Feature Importance
%s", display_title),
         x = "Predictor Variable", y = "Mean Decrease in Accuracy (Permutation)") +
    scale_fill_brewer(palette = "Set1") +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text.y      = element_text(size = 10, face = "bold"),
          legend.position  = "right")
  
  ggsave(file.path(RESULTS_DIR, sprintf("PFI_Color_30m_%s.png", tgt)), plot = p_color, width = 10, height = 6, dpi = 300)
  write_csv(imp_df %>% arrange(desc(Importance)), file.path(RESULTS_DIR, sprintf("PFI_Values_30m_%s.csv", tgt)))
}
