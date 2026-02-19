library(terra)
library(strucchange)
library(tictoc)

# --- CONFIGURACIÓN ---
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
  
  # 1. TEST DE ESTABILIDAD (Rápido)
  # Usamos un h ligeramente mayor para ganar velocidad en el escaneo
  qlr  <- strucchange::Fstats(ss ~ v_years, from = 0.05)
  test <- strucchange::sctest(qlr, type = "supF")
  
  f_stat  <- as.numeric(test$statistic)
  f_p_val <- as.numeric(test$p.value)
  
  # Inicializar resultados globales
  lm_total <- lm(ss ~ v_years)
  res_global <- c(
    n_breaks = 0,
    P_total  = as.numeric(coef(lm_total)[2]),
    Pv_total = summary(lm_total)$coefficients[2, 4],
    RSE      = sqrt(deviance(lm_total)/df.residual(lm_total)),
    T_mean_total = mean(ss),
    F_sup    = f_stat,
    P_val_sup = f_p_val
  )
  
  seg_data <- rep(NA, max_segs * 5)
  
  # 2. FILTRO DE SIGNIFICANCIA (La clave de la velocidad)
  # Solo buscamos los puntos exactos si el test dice que existen
  if (f_p_val < significance) {
    bp <- tryCatch({
      # Solo ejecutamos breakpoints si hay evidencia de cambio
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
    # Si no es significativo, el primer segmento es simplemente la serie completa
    seg_data[1:5] <- c(res_global["P_total"], res_global["Pv_total"], 1940, 1939 + n_total, res_global["T_mean_total"])
  }
  
  return(c(res_global, seg_data))
}

# --- PROCESAMIENTO ESPACIAL ---
r <- rast(path_nc)

# 1. LIMPIEZA DE PERIODOS INCOMPLETOS
# Supongamos que queremos solo años completos (1940-2025)
n_capas <- nlyr(r)
fechas <- seq(as.Date("1940-01-01"), by = "month", length.out = n_capas)
terra::time(r) <- fechas

# Filtramos: Solo nos quedamos con los meses de años menores a 2026
capas_validas <- which(as.numeric(format(fechas, "%Y")) <= 2025)
r <- r[[capas_validas]]

# 2. PROCESAMIENTO ANUAL
test_ext <- ext(10, 350, 35, 45) 
r <- crop(r, test_ext)
r_annual <- terra::tapp(r, "years", mean) - 273.15

# writeRaster(r_annual, "C:/A_TRABAJO/A_JORGE/ERA5/ERA5_anual_1940_2025.tif", overwrite = TRUE, datatype = "FLT4S")

tic()
res_rast <- terra::app(r_annual, fun = bp_pixel_analysis_fast, cores = 10)
toc()
final_df <- as.data.frame(res_rast, xy = TRUE)

# Nombres de columnas actualizados con el test F
base_names <- c("x", "y", "n_breaks", "P_total", "Pv_total", "RSE", "T_mean_total", "F_sup", "P_val_sup")
seg_names  <- c()
for(i in 1:max_segments_allowed) {
  seg_names <- c(seg_names, paste0(c("P_seg", "Pv_seg", "Start", "End", "Tmean_seg"), "_", i))
}
colnames(final_df) <- c(base_names, seg_names)
final_df

write.csv2(final_df, "C:/A_TRABAJO/A_JORGE/ERA5/RES_FINAL_CON_ESTADISTICO_F.csv")




# 1. Extraer la capa del número de rupturas (es la primera capa)
mapa_breaks <- res_rast[[1]]
names(mapa_breaks) <- "Numero_de_Rupturas"

# 2. Definir una paleta de colores discreta
# 0: gris/azul, 1: amarillo, 2: naranja, 3+: rojo
colores <- c("grey90", "#FFD700", "#FF8C00", "#FF4500", "#8B0000")

# 3. Dibujar el mapa
plot(mapa_breaks, 
     main = "Número de Puntos de Ruptura en la Tendencia (1940-2025)",
     col = colores,
     type = "classes", 
     levels = c(0, 1, 2, 3, 4),
     plg = list(title = "Rupturas", x = "bottomright"))

# Añadir fronteras para referencia (opcional)
# frontiers <- vect(algun_shapefile_de_paises)
# plot(frontiers, add = TRUE, border = "black")