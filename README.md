# TREM2_GFP_analysis
R scripts for automated quantification of spatial proximity between TREM2⁺ myeloid and GFP⁺ tumor cells in immunofluorescence images. Includes DAPI segmentation, intensity profiling, nearest-neighbour analysis, Ripley’s K/L statistics, and publication-ready outputs.
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

Results appear in /outputs/.
File	Description
MAP_seriesXX.png	Spatial map (DAPI, GFP⁺, TREM2⁺)
nearest_TREM2_to_GFP_ALL.xlsx	TREM2⁺→GFP⁺ distances (µm)
ripley_TREM2vGFP_ALL.csv	Ripley’s K/L spatial stats
recap_TREM2_vs_DAPI.xlsx	DAPI and TREM2⁺ counts per field
Citation

Author: Martin Pedard, University of Geneva
Contact: martin.pedard@unige.ch
]
Repository: https://github.com/pedarmar-git/TREM2_GFP_analysis

Year: 2025
