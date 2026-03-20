# TREM2–GFP Spatial Proximity Analysis

R scripts for quantifying spatial proximity between **TREM2⁺** and **GFP⁺** cells in immunofluorescence (IHF) images.

---

## Overview
This pipeline:
1. Segments DAPI nuclei (C=0).
2. Measures mean GFP (C=1 or C=2) and TREM2 (C=3) intensities per cell.
3. Classifies GFP⁺ and TREM2⁺ cells using adaptive thresholds.
4. Computes nearest-neighbour distances (TREM2⁺ → GFP⁺).
5. Calculates Ripley’s K/L spatial statistics.
6. Generates spatial maps (PNG), histograms, and Excel summaries.

---

## Requirements
- R ≥ 4.3  
- Packages: `EBImage`, `ggplot2`, `data.table`, `RANN`, `writexl`, `spatstat.geom`, `spatstat.explore`

---

## Usage
1. Place `.jpg` images in a folder:
(series 01) - C=0.jpg → DAPI
(series 01) - C=1.jpg → GFP
(series 01) - C=3.jpg → TREM2

2. Open `scripts/IHF_analysis_proximity_TREM2_GFP.R` in RStudio.  
3. Edit the path:
```r
img_dir <- "C:/Postdoctorat Lausanne/Manips/IF/Martin F480 TREM2/R IF TREM2/TREM2 IF"

4. Run:
source("scripts/IHF_analysis_proximity_TREM2_GFP.R")

5. Results appear in /outputs/.

Output files
For each processed field <id> (format: <base_root>__series<XX>):

File	Description
inventaire_champs.csv	List of all detected fields with per-channel file paths
MAP_<id>.png	Spatial map of cells colour-coded by status: grey = DAPI, green = GFP+, pink = TREM2+
cells_<id>.csv	Per-cell measurements for one field: coordinates, mean intensities, GFP/TREM2 status
nearest_TREM2_to_GFP_<id>.xlsx	TREM2+ → nearest GFP+ distances (double-positive self-pairing excluded) — per field
nearest_TREM2_to_GFP_<id>.csv	Same as above in CSV format
hist_TREM2_to_GFP_<id>.png	Histogram of TREM2+ → GFP+ distance distribution — per field
ripley_TREM2vGFP_<id>.csv	Ripley K/L cross-function values (TREM2 vs GFP) — per field
ripley_LminusR_<id>.png	Plot of L(r) − r from the cross Ripley function — per field
cells_ALL_bases_ALL_series.csv	Consolidated per-cell measurements across all fields
nearest_TREM2_to_GFP_ALL.xlsx	Consolidated TREM2→GFP distances across all fields
nearest_TREM2_to_GFP_ALL.csv	Same as above in CSV format
ripley_TREM2vGFP_ALL.csv	Consolidated Ripley K/L values across all fields
recap_TREM2_vs_DAPI.xlsx	Summary per field: n DAPI cells, n TREM2+, % TREM2+
recap_TREM2_vs_DAPI.csv	Same as above in CSV format

Author: Martin Pedard, University of Geneva
Contact: martin.pedard@unige.ch
]
Repository: https://github.com/pedarmar-git/TREM2_GFP_analysis

Year: 2025
