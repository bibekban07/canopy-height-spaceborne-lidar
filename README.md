# Cross-Platform GEDI and ICESat-2 Spaceborne LiDAR Integration and Multi-Sensor Satellite Data Fusion for Canopy Height Estimation in Fragmented Temperate Hardwood Forests

This repository provides open-source, reproducible R code and demonstration datasets accompanying the research manuscript submitted to *Remote Sensing*:

> **"Cross-Platform GEDI and ICESat-2 Spaceborne LiDAR Integration and Multi-Sensor Satellite Data Fusion for Canopy Height Estimation in Fragmented Temperate Hardwood Forests"**

---

## Repository Structure

```text
canopy-height-spaceborne-lidar/
|-- data/
|   |-- sample_training_pts.csv          # 5,000-row stratified demonstration sample
|   |-- sample_val15m.csv                # 15 m co-located sub-pixel holdout validation subset
|   |-- sample_val001deg.csv             # 0.001 deg stand-scale holdout validation subset
|   `-- sample_als_validation_pts.csv    # Independent airborne LiDAR (ALS) reference subset (1,647 points)
|-- scripts/
|   |-- 01_semivariogram_spatial_blocks.R           # Empirical variogram modeling justifying spatial block sizes
|   |-- 02_master_fusion_and_four_tier_validation.R # 4-tier validation pipeline (Random, Spatial, 15m, 0.001 deg)
|   |-- 03_permutation_feature_importance.R         # Permutation Feature Importance (PFI) across 5 domains
|   `-- 04_als_independent_validation.R             # Independent ALS point validation and plotting
|-- results/                             # Output tables, validation matrices, and publication figures
|-- .gitignore
`-- README.md
```

---

## Demonstration Data Note
The CSV files in data/ are representative demonstration subsets sampled across Southern Illinois. They allow immediate execution and functional verification of the complete modeling and validation pipeline on standard workstations without requiring multi-gigabyte spatial raster downloads.

Regional wall-to-wall multi-spectral satellite stacks (Landsat 8, Sentinel-1 SAR, Sentinel-2 Red-Edge, Topography) and continuous Airborne Laser Scanning (ALS) canopy height models are permanently archived under institutional repositories upon final publication.

---

## Reproduction Instructions
Open R or RStudio with the working directory set to the project root, then execute any script:
```R
source("scripts/01_semivariogram_spatial_blocks.R")
source("scripts/02_master_fusion_and_four_tier_validation.R")
source("scripts/03_permutation_feature_importance.R")
source("scripts/04_als_independent_validation.R")
```
