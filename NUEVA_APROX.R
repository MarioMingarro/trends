library(terra)
library(strucchange)
library(sandwich)
library(tictoc)

# --- 1. CARGA Y PREPARACIÓN DE DATOS ---
path_nc <- "C:/A_TRABAJO/A_JORGE/ERA5/ERA5_DATA/era5_1940_2025.nc"
r <- rast(path_nc)

# Configuración temporal
n_capas <- nlyr(r)
fechas <- seq(as.Date("1940-01-01"), by = "month", length.out = n_capas)
terra::time(r) <- fechas
capas_validas <- which(as.numeric(format(fechas, "%Y")) <= 2025)
r <- r[[capas_validas]]

# Área de estudio y paso a anual (Celsius)
test_ext <- ext(0, 10, 35, 45) 
r_annual <- terra::crop(r, test_ext)
r_annual <- terra::tapp(r_annual, "years", mean) - 273.15
r_annual <- rast("C:/A_TRABAJO/A_JORGE/ERA5/ERA5_anual_1940_2025.tif")
# --- 2. FUNCIÓN DE ANÁLISIS MEJORADA (con Aceleración) ---
bp_pixel_analysis_final <- function(x) {
  n_total <- length(x)
  v_years <- 1:n_total
  ss <- as.numeric(x)
  
  # 18 columnas de salida (incluyendo 2 de aceleración)
  out <- rep(NA, 18)
  
  if (any(is.na(x)) || length(unique(na.omit(x))) < 30) return(out)
  
  # 1. TEST GLOBAL HAC
  qlr <- strucchange::Fstats(ss ~ v_years, from = 0.15)
  test_hac <- try(strucchange::sctest(qlr, type = "supF", vcov = sandwich::vcovHAC), silent = TRUE)
  if(inherits(test_hac, "try-error")) return(out)
  
  out[1] <- as.numeric(test_hac$statistic)
  out[2] <- as.numeric(test_hac$p.value)
  
  # 2. DETECCIÓN Y TENDENCIAS SEGMENTADAS
  if (!is.na(out[2]) && out[2] < 0.05) {
    bp <- try(strucchange::breakpoints(ss ~ v_years, h = 0.15, breaks = 2), silent = TRUE)
    
    if (!inherits(bp, "try-error") && !all(is.na(bp$breakpoints))) {
      actual_breaks <- sort(bp$breakpoints)
      n_breaks <- length(actual_breaks)
      
      out[3] <- 1939 + actual_breaks[1] # Year_1
      if(n_breaks > 1) out[4] <- 1939 + actual_breaks[2] # Year_2
      
      z1 <- pmax(0, v_years - actual_breaks[1])
      
      if(n_breaks == 1) {
        fit <- lm(ss ~ v_years + z1)
        s <- summary(fit)$coefficients
        if(nrow(s) >= 3){
          out[5] <- s[2,1]           # Tr1
          out[6] <- s[2,1] + s[3,1]  # Tr2
          out[8] <- s[2,4]           # Pv_Tr1
          out[9] <- s[3,4]           # Pv_C1 (Significancia del cambio)
          # ACELERACIÓN 1: Diferencia de pendientes
          out[17] <- s[3,1]          # Acc_1 (Tr2 - Tr1)
        }
      } else {
        z2 <- pmax(0, v_years - actual_breaks[2])
        fit <- lm(ss ~ v_years + z1 + z2)
        s <- summary(fit)$coefficients
        if(nrow(s) >= 4){
          out[5] <- s[2,1]                    # Tr1
          out[6] <- s[2,1] + s[3,1]           # Tr2
          out[7] <- s[2,1] + s[3,1] + s[4,1]  # Tr3
          out[8] <- s[2,4]   # Pv_Tr1
          out[9] <- s[3,4]   # Pv_C1
          out[10] <- s[4,4]  # Pv_C2
          # ACELERACIONES:
          out[17] <- s[3,1]  # Acc_1 (Tr2 - Tr1)
          out[18] <- s[4,1]  # Acc_2 (Tr3 - Tr2)
        }
      }
      
      # Medias por segmento
      lims <- c(0, actual_breaks, n_total)
      for(i in 1:(length(lims)-1)) {
        if(i <= 3) out[10+i] <- mean(ss[(lims[i]+1):lims[i+1]], na.rm=TRUE)
      }
    }
  }
  
  # 3. TENDENCIA GLOBAL
  fit_glob <- try(lm(ss ~ v_years), silent = TRUE)
  if(!inherits(fit_glob, "try-error")) {
    s_glob <- summary(fit_glob)$coefficients
    out[14] <- s_glob[2,1] # Global_Tr
    out[15] <- s_glob[2,4] # Global_Pv
    out[16] <- mean(ss, na.rm=TRUE) # Global_Tm
  }
  return(out)
}

# --- 3. EJECUCIÓN ---
tic("Análisis Estructural con Aceleración")
res_rast <- terra::app(r_annual, fun = bp_pixel_analysis_final, cores = 10)
toc()

# --- 4. EXTRACCIÓN Y LIMPIEZA ---
final_df <- as.data.frame(res_rast, xy = TRUE)
colnames(final_df) <- c(
  "x", "y", "F_HAC", "Pv_HAC", "Year_1", "Year_2",
  "Tr_1", "Tr_2", "Tr_3", "Pv_Tr1", "Pv_C1", "Pv_C2",
  "Tm_1", "Tm_2", "Tm_3", "Global_Tr", "Global_Pv", "Global_Tm",
  "Acc_1", "Acc_2"
)
write.csv2(final_df, "C:/A_TRABAJO/A_JORGE/ERA5/RES_FINAL_NUEVA_2026.csv")
