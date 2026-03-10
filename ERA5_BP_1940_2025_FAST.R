library(terra)
library(strucchange)
library(tictoc)


path_nc <- "C:/A_TRABAJO/A_JORGE/ERA5/ERA5_DATA/era5_1940_2025.nc"
max_segments_allowed <- 5 


bp_pixel_analysis_fast <- function(x) {
  max_segs <- 5 
  significance <- 0.01
  n_total <- length(x)

  if (any(is.na(x)) || all(x == 0) || length(unique(na.omit(x))) < 10) {
    return(rep(NA, 7 + (max_segs * 5)))
  }
  
  v_years <- 1:n_total
  ss <- ts(x, frequency = 1)
  
  # 1. TEST DE ESTABILIDAD (SupF Test)
  # h=0.05 permite capturar segmentos de hasta ~4 años
  qlr  <- strucchange::Fstats(ss ~ v_years, from = 0.05)
  test <- strucchange::sctest(qlr, type = "supF")
  
  f_stat  <- as.numeric(test$statistic)
  f_p_val <- as.numeric(test$p.value)
  

  lm_total <- lm(ss ~ v_years)
  res_global <- c(
    n_breaks     = 0,
    P_total      = as.numeric(coef(lm_total)[2]),
    Pv_total     = summary(lm_total)$coefficients[2, 4],
    RSE          = sqrt(deviance(lm_total)/df.residual(lm_total)),
    T_mean_total = mean(ss),
    F_sup        = f_stat,
    P_val_sup    = f_p_val
  )

  seg_data <- rep(NA, max_segs * 5)
  
  # 2. LÓGICA DE RUPTURAS
  # Si el test F detecta inestabilidad, buscamos dónde están los quiebres
  if (f_p_val < significance) {
    bp <- tryCatch({
      # Dejamos que el algoritmo encuentre el número óptimo de rupturas
      strucchange::breakpoints(ss ~ v_years, h = 0.05) 
    }, error = function(e) return(NULL))
    
    if (!is.null(bp) && !is.na(bp$breakpoints[1])) {
      break_pts <- bp$breakpoints
      res_global["n_breaks"] <- length(break_pts)
      
      pts <- c(1, break_pts, n_total)
      n_segments <- length(pts) - 1
    for (j in 1:min(n_segments, max_segs)) {
        idx_start <- if (j == 1) pts[j] else pts[j] + 1
        idx_end   <- pts[j+1]
        
        if (idx_start < idx_end) {
          y_seg   <- v_years[idx_start:idx_end]
          val_seg <- ss[idx_start:idx_end]
          t_med_seg <- mean(val_seg)
          
          p_seg <- NA; pv_seg <- NA
          if (length(val_seg) > 3) {
            lm_seg <- lm(val_seg ~ y_seg)
            p_seg  <- coef(lm_seg)[2]
            pv_seg <- summary(lm_seg)$coefficients[2, 4]
          }
          
          pos <- (j - 1) * 5 + 1
          seg_data[pos:(pos+4)] <- c(p_seg, pv_seg, 1939 + idx_start, 1939 + idx_end, t_med_seg)
        }
      }
    }
  } else {
    seg_data[1:5] <- c(res_global["P_total"], res_global["Pv_total"], 1940, 1939 + n_total, res_global["T_mean_total"])
  }
  
  return(c(res_global, seg_data))
}

# --- PROCESAMIENTO ESPACIAL ---
r <- rast(path_nc)


n_capas <- nlyr(r)
fechas <- seq(as.Date("1940-01-01"), by = "month", length.out = n_capas)
terra::time(r) <- fechas


capas_validas <- which(as.numeric(format(fechas, "%Y")) <= 2025)
r <- r[[capas_validas]]

# EXTEND TEST
test_ext <- ext(0, 10, 35, 45) 
r <- crop(r, test_ext)

# Media anual y conversión de Kelvin a Celsius
r_annual <- terra::tapp(r, "years", mean) - 273.15


tic("Procesamiento de tendencias")
res_rast <- terra::app(r_annual, fun = bp_pixel_analysis_fast, cores = 10)
toc()

# --- EXTRACCIÓN Y ETIQUETADO ---
final_df <- as.data.frame(res_rast, xy = TRUE)

base_names <- c("n_breaks", "P_total", "Pv_total", "RSE", "T_mean_total", "F_sup", "P_val_sup")

seg_vars <- c("P_seg", "Pv_seg", "Start", "End", "Tmean_seg")
seg_names <- expand.grid(seg_vars, 1:max_segments_allowed)
seg_names <- paste0(seg_names$Var1, "_", seg_names$Var2)

colnames(final_df) <- c("x", "y", base_names, seg_names)

# --- RESULTADOS ---
head(final_df)
