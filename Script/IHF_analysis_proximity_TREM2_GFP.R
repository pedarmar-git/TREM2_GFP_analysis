# ============ IHF — TREM2→GFP nearest (no self), cartes, Excel, histogrammes, Ripley K/L ============
suppressPackageStartupMessages({
  library(EBImage); library(ggplot2); library(data.table); library(stringr); library(grid)
})
suppressWarnings({
  if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
  if (!requireNamespace("RANN", quietly = TRUE))     install.packages("RANN")
  if (!requireNamespace("spatstat.geom", quietly = TRUE))    install.packages("spatstat.geom")
  if (!requireNamespace("spatstat.explore", quietly = TRUE)) install.packages("spatstat.explore")
})
library(writexl); library(RANN)
library(spatstat.geom); library(spatstat.explore)

options(stringsAsFactors = FALSE)

# -------- PARAMÈTRES --------
img_dir <- "C:\Users\pedard\Desktop\Github deposit"
out_dir <- file.path(img_dir, "outputs"); if (!dir.exists(out_dir)) dir.create(out_dir, TRUE)

# Vitesse & segmentation
SPEED_MODE    <- TRUE
USE_WATERSHED <- FALSE
TARGET_N      <- 1200L
TARGET_RANGE  <- c(700L, 1800L)

# Mesure (anneaux distincts)
MEASURE_RING_GFP   <- TRUE
MEASURE_RING_TREM2 <- TRUE
RING_RADIUS_GFP    <- 4
RING_RADIUS_TREM2  <- 2

# Soustraction de fond
GFP_BG_RADIUS   <- 31
TREM2_BG_RADIUS <- 41

# Seuils adaptatifs
GFP_PREF    <- "otsu"
GFP_QFALL   <- 0.40
TREM2_PREF  <- "otsu"
TREM2_QFALL <- 0.85

# Distance / proximité
PIXEL_SIZE_UM  <- 0.5              # µm / pixel
PROX_THRESH_UM <- c(10, 20, 50)

# Rendu cartes
dot_radius  <- 2.6
dot_stroke  <- 0.35
dot_scale   <- 1.25
dapi_fill   <- "#8A8A8A"
dapi_alpha  <- 0.35
gfp_alpha   <- 0.80
trem2_alpha <- 0.95

# Histogramme distances
HIST_BINWIDTH_UM <- 2
HIST_XMAX_UM     <- 60
HIST_BASESIZE    <- 12

# Ripley K/L cross (TREM2 vs GFP)
RIPLEY_R_MAX_UM  <- 100           # rayon max
RIPLEY_R_STEP_UM <- 2             # pas
RIPLEY_CORR      <- "isotropic"   # correction de bord

# -------- UTILS --------
logm <- function(...) { message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...)); flush.console() }
to_gray <- function(x){
  if (length(dim(x))==3L) x <- EBImage::channel(x,"gray")
  EBImage::colorMode(x) <- EBImage::Grayscale
  EBImage::normalize(x)
}
bg_subtract <- function(img, radius = 31) {
  img <- to_gray(img)
  bg  <- EBImage::opening(img, EBImage::makeBrush(radius, "disc"))
  EBImage::normalize(pmax(img - bg, 0))
}
pick_threshold <- function(v, prefer="otsu", q_fallback=0.45){
  v <- as.numeric(v); v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  if (prefer=="otsu") {
    br <- max(64L, round(sqrt(length(v)))); h <- hist(v, breaks=br, plot=FALSE)
    thr <- tryCatch({
      mids <- h$mids; counts <- h$counts; p <- counts/sum(counts)
      omega <- cumsum(p); mu <- cumsum(p*mids); mu_t <- mu[length(mu)]
      sigma_b2 <- (mu_t*omega - mu)^2 / (omega*(1-omega) + 1e-12)
      mids[which.max(sigma_b2)]
    }, error=function(e) NA_real_)
    if (is.finite(thr)) return(thr)
  }
  suppressWarnings(as.numeric(stats::quantile(v, probs=q_fallback, na.rm=TRUE)))
}

# -------- SEGMENTATION DAPI (FAST) --------
segment_dapi_fast <- function(dapi){
  d <- to_gray(dapi)
  if (SPEED_MODE) { sigmas <- 6; wins <- 23; offs <- 0.00 }
  else            { sigmas <- c(4,6,8); wins <- c(19,23,27); offs <- c(-0.02,0.00,0.02) }
  
  best <- list(score=Inf, lab=NULL)
  for(sg in sigmas) for(w in wins) for(off in offs){
    x <- EBImage::normalize(d - EBImage::gblur(d, sigma=sg)); x[x<0] <- 0; x <- EBImage::medianFilter(x, 3)
    bw <- EBImage::thresh(x, w, w, off) | (x > EBImage::otsu(x))
    bw <- EBImage::opening(bw, EBImage::makeBrush(3,"disc"))
    bw <- EBImage::closing(bw, EBImage::makeBrush(3,"disc"))
    bw <- EBImage::fillHull(bw)
    if (USE_WATERSHED) {
      dist <- EBImage::distmap(bw)
      seeds <- EBImage::watershed(1 - EBImage::normalize(dist), tolerance = 1, ext = 1)
      bw <- EBImage::propagate(x, seeds, mask = bw) > 0
    }
    lab <- EBImage::bwlabel(bw)
    if (max(lab)==0) next
    shp <- EBImage::computeFeatures.shape(lab)
    keep <- which(shp[,"s.area"] >= 10 & shp[,"s.area"] <= 8000)
    if (!length(keep)) next
    lab <- EBImage::rmObjects(lab, setdiff(seq_len(max(lab)), keep))
    n <- as.integer(max(lab))
    score <- abs(n - TARGET_N)
    if (score < best$score) best <- list(score=score, lab=EBImage::bwlabel(lab > 0))
    if (n >= TARGET_RANGE[1] && n <= TARGET_RANGE[2]) return(best$lab)
    if (SPEED_MODE) break
  }
  if (is.null(best$lab)) EBImage::bwlabel(0*d) else best$lab
}

# -------- MESURES (anneaux distincts GFP/TREM2) --------
measure_cells_dual <- function(lab, gfp, trem2,
                               ring_gfp_px = 0, use_ring_gfp = FALSE,
                               ring_trem2_px = 0, use_ring_trem2 = FALSE){
  if (is.null(lab) || max(lab) == 0) return(data.table())
  
  lab_g <- lab
  if (use_ring_gfp && ring_gfp_px > 0) {
    mask_nuc <- lab > 0
    mask_exp <- EBImage::dilate(mask_nuc, EBImage::makeBrush(2*ring_gfp_px + 1, "disc"))
    lab_g    <- EBImage::propagate(1 - EBImage::normalize(EBImage::distmap(mask_nuc)), seeds = lab, mask = mask_exp)
  }
  lab_r <- lab
  if (use_ring_trem2 && ring_trem2_px > 0) {
    mask_nuc <- lab > 0
    mask_exp <- EBImage::dilate(mask_nuc, EBImage::makeBrush(2*ring_trem2_px + 1, "disc"))
    lab_r    <- EBImage::propagate(1 - EBImage::normalize(EBImage::distmap(mask_nuc)), seeds = lab, mask = mask_exp)
  }
  
  shp <- EBImage::computeFeatures.shape(lab);  if (is.null(dim(shp)) || nrow(shp)==0) return(data.table())
  mom <- EBImage::computeFeatures.moment(lab); if (is.null(dim(mom)) || nrow(mom)==0) return(data.table())
  fg  <- EBImage::computeFeatures.basic(lab_g, gfp);    if (is.null(dim(fg)))  return(data.table())
  fr  <- EBImage::computeFeatures.basic(lab_r, trem2);  if (is.null(dim(fr)))  return(data.table())
  
  shp <- as.data.frame(shp); shp$id <- as.integer(rownames(shp)); shp <- shp[, c("id","s.area")]
  mom <- as.data.frame(mom); mom$id <- as.integer(rownames(mom)); mom <- mom[, c("id","m.cx","m.cy")]
  fg  <- as.data.frame(fg);  fg$id  <- as.integer(rownames(fg));  names(fg)[names(fg)=="b.mean"] <- "mean_gfp";   fg <- fg[, c("id","mean_gfp")]
  fr  <- as.data.frame(fr);  fr$id  <- as.integer(rownames(fr));  names(fr)[names(fr)=="b.mean"] <- "mean_trem2"; fr <- fr[, c("id","mean_trem2")]
  
  dt <- Reduce(function(a,b) merge(a,b, by="id", all=FALSE), list(shp, mom, fg, fr))
  if (!nrow(dt)) return(data.table())
  setDT(dt); setnames(dt, c("m.cx","m.cy","s.area"), c("x","y","area_px"))
  dt[, .(id, x, y, area_px, mean_gfp, mean_trem2)]
}

# -------- Proximité / mapping (exclusion du "self" si double-positif) --------
nearest_map_excl_self <- function(dt, px_um = 0.5, from_flag = "rpos", to_flag = "gpos"){
  A <- dt[get(from_flag) == TRUE, .(id, x, y)]
  B <- dt[get(to_flag)   == TRUE, .(id, x, y)]
  if (!nrow(A) || !nrow(B)) {
    return(data.table(from_id=integer(0), from_x=numeric(0), from_y=numeric(0),
                      to_id=integer(0),   to_x=numeric(0),   to_y=numeric(0),
                      dist_um=numeric(0)))
  }
  # exclure auto-correspondance si l'ID est commun (double-positives)
  inter_ids <- intersect(A$id, B$id)
  if (length(inter_ids)) B <- B[!id %in% inter_ids]
  if (!nrow(B)) {
    return(data.table(from_id=A$id, from_x=A$x, from_y=A$y,
                      to_id=NA_integer_, to_x=NA_real_, to_y=NA_real_,
                      dist_um=NA_real_))
  }
  am <- as.matrix(A[, .(x, y)])
  bm <- as.matrix(B[, .(x, y)])
  nn <- RANN::nn2(data = bm, query = am, k = 1)
  idx <- as.integer(nn$nn.idx)
  dx  <- am[,1] - bm[idx,1]; dy <- am[,2] - bm[idx,2]
  d_um <- sqrt(dx*dx + dy*dy) * px_um
  data.table(
    from_id = A$id, from_x = A$x, from_y = A$y,
    to_id   = B$id[idx], to_x = B$x[idx], to_y = B$y[idx],
    dist_um = d_um
  )
}

# -------- Ripley K/L (cross TREM2 vs GFP) --------
ripley_cross <- function(dt, px_um = 0.5, r_max_um = 100, r_step_um = 2, corr = "isotropic", w_px=NULL, h_px=NULL){
  TREM2 <- dt[rpos == TRUE, .(x_um = x*px_um, y_um = y*px_um)]
  GFP   <- dt[gpos == TRUE, .(x_um = x*px_um, y_um = y*px_um)]
  if (!nrow(TREM2) || !nrow(GFP)) return(data.table())
  if (is.null(w_px) || is.null(h_px)) {
    w_px <- max(dt$x, na.rm=TRUE); h_px <- max(dt$y, na.rm=TRUE)
  }
  W <- spatstat.geom::owin(xrange = c(0, w_px*px_um), yrange = c(0, h_px*px_um))
  X <- rbind(TREM2[, .(x_um, y_um, mark="TREM2")], GFP[, .(x_um, y_um, mark="GFP")])
  p <- spatstat.geom::ppp(x = X$x_um, y = X$y_um, window = W,
                          marks = factor(X$mark, levels=c("TREM2","GFP")),
                          unitname = c("micrometre","micrometres"))
  rseq <- seq(0, r_max_um, by = r_step_um)
  kcr  <- spatstat.explore::Kcross(p, i="TREM2", j="GFP", correction = corr, r = rseq)
  corr_col <- switch(corr, "isotropic"="iso", "border"="bord", "translate"="trans", "iso")
  K <- as.data.table(kcr)[, .(r = r, K = get(corr_col))]
  K[, L := sqrt(K/pi)]
  K[, LminusR := L - r]
  setnames(K, c("r","K","L","LminusR"), c("r_um","Kcross","Lcross","LminusR"))
  K
}

# -------- PARSING FICHIERS --------
stopifnot(dir.exists(img_dir))
cat("\n=== Vérification du dossier ===\n")
files_test <- list.files(img_dir, pattern="\\.(jpg|jpeg)$", full.names=TRUE, recursive=TRUE, ignore.case=TRUE)
if (!length(files_test)) stop("Aucun .jpg/.jpeg trouvé — corrigez 'img_dir' ou placez vos images.")
cat(sprintf("✅ %d fichiers image trouvés (aperçu):\n", length(files_test)))
print(utils::head(basename(files_test), 5))
cat("===============================\n\n"); flush.console()

files <- list.files(img_dir, pattern="\\.(jpg|jpeg)$", full.names=TRUE, ignore.case=TRUE, recursive=TRUE)
dtf <- data.table(path=files, file=basename(files))
m <- stringr::str_match(dtf$file, "^(.*)\\s\\(\\s*series\\s*(\\d+)\\s*\\)\\s-\\sC\\s*=\\s*(\\d+)\\.(?:jpg|jpeg)$")
dtf[, base_root := tolower(trimws(m[,2]))]
dtf[, series_id := suppressWarnings(as.integer(m[,3]))]
dtf[, chan      := m[,4]]
dtf <- dtf[!is.na(base_root) & !is.na(series_id) & !is.na(chan)]
if (!nrow(dtf)) stop("Parsing impossible : attendu '... (series XX) - C=Y.jpg'")

grp <- dtf[, .(
  pC0 = path[chan=="0"][1],
  pC1 = path[chan=="1"][1],
  pC2 = path[chan=="2"][1],
  pC3 = path[chan=="3"][1]
), by=.(base_root, series_id)]
grp <- grp[!is.na(pC0) & !is.na(pC3) & (!is.na(pC1) | !is.na(pC2))]
if (!nrow(grp)) stop("Pas de triplet complet (C0 DAPI, C3 TREM2, C1/C2 GFP).")
setorder(grp, base_root, series_id)
fwrite(grp, file.path(out_dir, "inventaire_champs.csv"))
logm("Champs détectés: %d (bases=%d, séries uniques=%d)", nrow(grp), uniqueN(grp$base_root), uniqueN(grp$series_id))

# -------- TRAITEMENT --------
stash_cells <- list()
stash_map   <- list()
stash_ripley<- list()
summary_rows <- list()

for (i in seq_len(nrow(grp))) {
  r <- grp[i]
  base <- r$base_root; sid <- r$series_id
  id_lab <- sprintf("%s__series%02d", gsub("[^a-z0-9]+","_", base), as.integer(sid))
  logm("=== Champ %d/%d: base='%s' | series=%02d ===", i, nrow(grp), base, as.integer(sid))
  
  p_dapi  <- r$pC0; p_trem2 <- r$pC3; p_gfp <- if (!is.na(r$pC2)) r$pC2 else r$pC1
  dapi  <- tryCatch(EBImage::readImage(p_dapi),  error=function(e){ logm("  ❌ read DAPI: %s",  e$message); NULL })
  gfp   <- tryCatch(EBImage::readImage(p_gfp),   error=function(e){ logm("  ❌ read GFP: %s",   e$message); NULL })
  trem2 <- tryCatch(EBImage::readImage(p_trem2), error=function(e){ logm("  ❌ read TREM2: %s", e$message); NULL })
  if (any(sapply(list(dapi,gfp,trem2), is.null))) { logm("  → champ ignoré (lecture échouée)."); next }
  
  lab <- segment_dapi_fast(dapi)
  if (max(lab)==0) { logm("  ⚠️ Segmentation vide → ignoré."); next }
  
  dapi_g  <- to_gray(dapi); gfp_g <- to_gray(gfp); trem2_g <- to_gray(trem2)
  common  <- dim(dapi_g)[1:2]
  gfp_g   <- EBImage::resize(gfp_g,   w=common[2], h=common[1])
  trem2_g <- EBImage::resize(trem2_g, w=common[2], h=common[1])
  
  gfp_bs   <- bg_subtract(gfp_g,   GFP_BG_RADIUS)
  trem2_bs <- bg_subtract(trem2_g, TREM2_BG_RADIUS)
  
  dtc <- measure_cells_dual(
    lab, gfp_bs, trem2_bs,
    ring_gfp_px   = RING_RADIUS_GFP,   use_ring_gfp   = MEASURE_RING_GFP,
    ring_trem2_px = RING_RADIUS_TREM2, use_ring_trem2 = MEASURE_RING_TREM2
  )
  if (!nrow(dtc)) { logm("  ⚠️ Aucune cellule mesurée → ignoré."); next }
  
  thr_gfp <- min(
    pick_threshold(dtc$mean_gfp,   prefer=GFP_PREF,   q_fallback=GFP_QFALL),
    suppressWarnings(as.numeric(stats::quantile(dtc$mean_gfp,   probs=0.45, na.rm=TRUE)))
  )
  thr_trem2 <- max(
    pick_threshold(dtc$mean_trem2, prefer=TREM2_PREF, q_fallback=TREM2_QFALL),
    suppressWarnings(as.numeric(stats::quantile(dtc$mean_trem2, probs=0.90, na.rm=TRUE)))
  )
  
  dtc[, gpos := (mean_gfp   >= thr_gfp)]
  dtc[, rpos := (mean_trem2 >= thr_trem2)]
  
  # -------- CARTE SPATIALE --------
  dot_radius2 <- dot_radius * dot_scale
  pmap <- ggplot() +
    geom_point(data = dtc, aes(x = x, y = y),
               shape = 21, size = dot_radius2, stroke = 0,
               fill  = dapi_fill, alpha = dapi_alpha) +
    geom_point(data = dtc[gpos == TRUE & rpos == FALSE], aes(x = x, y = y),
               shape = 21, size = dot_radius2, stroke = dot_stroke,
               fill  = "#2ECC71", color = "#1B8F4E", alpha = gfp_alpha) +
    geom_point(data = dtc[rpos == TRUE], aes(x = x, y = y),
               shape = 21, size = dot_radius2, stroke = dot_stroke,
               fill  = "#D62780", color = "#7A1C4C", alpha = trem2_alpha) +
    coord_fixed() + scale_y_reverse() +
    theme_void(base_size = 12) +
    theme(plot.background  = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA))
  ggsave(file.path(out_dir, paste0("MAP_", id_lab, ".png")), pmap, width = 7, height = 7, dpi = 300)
  
  # -------- DISTANCES TREM2+ → GFP+ (exclusion de soi) --------
  map_r2g <- nearest_map_excl_self(dtc, px_um = PIXEL_SIZE_UM, from_flag = "rpos", to_flag = "gpos")
  if (nrow(map_r2g)) map_r2g[, `:=`(base_root = base, series_id = sprintf("%02d", as.integer(sid)))]
  
  writexl::write_xlsx(map_r2g, file.path(out_dir, paste0("nearest_TREM2_to_GFP_", id_lab, ".xlsx")))
  data.table::fwrite(map_r2g,  file.path(out_dir, paste0("nearest_TREM2_to_GFP_", id_lab, ".csv")))
  logm("  ✅ Excel distances (TREM2→GFP, no self) écrit.")
  
  # -------- HISTOGRAMME DISTANCES --------
  if (nrow(map_r2g)) {
    ph <- ggplot(map_r2g[is.finite(dist_um)], aes(x = dist_um)) +
      geom_histogram(binwidth = HIST_BINWIDTH_UM, boundary = 0, closed = "left") +
      coord_cartesian(xlim = c(0, HIST_XMAX_UM)) +
      theme_minimal(base_size = HIST_BASESIZE) +
      labs(x = "Distance TREM2+ → GFP+ (µm)", y = "Nb de cellules",
           title = paste0("Distribution des distances — ", gsub("_"," ", id_lab)))
    ggsave(file.path(out_dir, paste0("hist_TREM2_to_GFP_", id_lab, ".png")), ph, width = 6, height = 4, dpi = 300)
  }
  
  # -------- RIPLEY K/L (TREM2 vs GFP) --------
  rip <- ripley_cross(dtc, px_um = PIXEL_SIZE_UM,
                      r_max_um = RIPLEY_R_MAX_UM, r_step_um = RIPLEY_R_STEP_UM,
                      corr = RIPLEY_CORR,
                      w_px = common[2], h_px = common[1])
  if (nrow(rip)) {
    rip[, `:=`(base_root = base, series_id = sprintf("%02d", as.integer(sid)))]
    data.table::fwrite(rip, file.path(out_dir, paste0("ripley_TREM2vGFP_", id_lab, ".csv")))
    pl <- ggplot(rip, aes(x = r_um, y = LminusR)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_line() +
      theme_minimal(base_size = 12) +
      labs(x = "r (µm)", y = "L(r) - r",
           title = paste0("Ripley Lcross (TREM2→GFP) — ", gsub("_"," ", id_lab)))
    ggsave(file.path(out_dir, paste0("ripley_LminusR_", id_lab, ".png")), pl, width = 6, height = 4, dpi = 300)
  }
  
  # -------- Stash globaux --------
  data.table::fwrite(dtc[, .(id, x, y, area_px, mean_gfp, mean_trem2, gpos, rpos)],
                     file.path(out_dir, paste0("cells_", id_lab, ".csv")))
  stash_cells[[length(stash_cells)+1]] <- copy(dtc)[, `:=`(base_root = base, series_id = sprintf("%02d", as.integer(sid)))]
  stash_map[[length(stash_map)+1]]     <- map_r2g
  if (nrow(rip)) stash_ripley[[length(stash_ripley)+1]] <- rip
  
  # Récap DAPI/TREM2 par champ
  summary_rows[[length(summary_rows)+1]] <- data.table(
    base_root  = base,
    series_id  = sprintf("%02d", as.integer(sid)),
    n_DAPI     = nrow(dtc),
    n_TREM2pos = sum(dtc$rpos, na.rm = TRUE),
    pct_TREM2  = round(100*sum(dtc$rpos, na.rm = TRUE)/max(1, nrow(dtc)), 1)
  )
}

# -------- EXPORTS GLOBAUX --------
if (length(stash_cells)) {
  all_cells <- data.table::rbindlist(stash_cells, use.names=TRUE, fill=TRUE)
  data.table::fwrite(all_cells, file.path(out_dir, "cells_ALL_bases_ALL_series.csv"))
  logm("✅ CSV global cellules écrit.")
}
if (length(stash_map)) {
  all_map <- data.table::rbindlist(stash_map, use.names=TRUE, fill=TRUE)
  data.table::setorder(all_map, base_root, series_id, from_id)
  writexl::write_xlsx(all_map, file.path(out_dir, "nearest_TREM2_to_GFP_ALL.xlsx"))
  data.table::fwrite(all_map,        file.path(out_dir, "nearest_TREM2_to_GFP_ALL.csv"))
  logm("✅ Excel/CSV global distances (TREM2→GFP, no self) écrits.")
}
if (length(stash_ripley)) {
  all_rip <- data.table::rbindlist(stash_ripley, use.names=TRUE, fill=TRUE)
  data.table::fwrite(all_rip, file.path(out_dir, "ripley_TREM2vGFP_ALL.csv"))
  logm("✅ CSV global Ripley K/L écrit.")
}
if (length(summary_rows)) {
  recap <- data.table::rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
  data.table::setorder(recap, base_root, series_id)
  writexl::write_xlsx(recap, file.path(out_dir, "recap_TREM2_vs_DAPI.xlsx"))
  data.table::fwrite(recap,      file.path(out_dir, "recap_TREM2_vs_DAPI.csv"))
  logm("✅ Récap DAPI/TREM2 écrit.")
}

logm("✅ Terminé. Sorties dans: %s", out_dir)
