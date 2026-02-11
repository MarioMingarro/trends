library(raster)
library(terra)
library(strucchange)
library(doParallel)
library(foreach)
library(tictoc)

# --- 1. Definir la función ANTES de usarla ---
bp_analysis <- function(i, data, start_year, end_year) {
  # Obtener los datos de la serie temporal (quitando columnas lon/lat)
  ss_vals <- as.numeric(data[i, -c(1, 2)])
  
  # Manejo de NAs y longitud de serie
  if(any(is.na(ss_vals))) {
    # Retornar una fila vacía con la misma estructura si hay NAs
    return(data.frame(matrix(NA, nrow=1, ncol=48))) # Ajustar ncol a los resultados
  }
  
  ss <- ts(ss_vals, start = start_year, end = end_year, frequency = 1)
  
  # Detección de rupturas
  bp <- strucchange::breakpoints(ss ~ 1)
  breakdates <- strucchange::breakdates(bp)
  
  # Inicializar resultados
  results <- list(
    year_break_1 = 999, year_break_2 = 999, year_break_3 = 999, year_break_4 = 999, year_break_5 = 999,                
    n_breaks = length(breakdates), F_sup = NA, p_value = NA,
    P_total = NA, Pv_total = NA, RSE_total = NA, tmed_total = NA
  )
  
  # Inicializar segmentos pre/post hasta 6 segmentos posibles (5 rupturas)
  for(k in 1:6) {
    results[[paste0("P_pre_", k)]] <- NA
    results[[paste0("Pv_pre_", k)]] <- NA
    results[[paste0("RSE_pre_", k)]] <- NA
    results[[paste0("tmed_pre_", k)]] <- NA
    results[[paste0("P_post_", k)]] <- NA
    results[[paste0("Pv_post_", k)]] <- NA
    results[[paste0("RSE_post_", k)]] <- NA
    results[[paste0("tmed_post_", k)]] <- NA
  }
  
  # Análisis estadístico de rupturas
  qlr <- strucchange::Fstats(ss ~ 1, from = 0.05) # 'from' ajustado para series cortas
  test <- strucchange::sctest(qlr, type = "supF")
  results$F_sup <- test$statistic
  results$p_value <- test$p.value
  
  # Análisis total
  years_seq <- seq(start_year, end_year, 1)
  lm_total <- lm(ss ~ years_seq)
  results$P_total <- round(coef(lm_total)[2], 4)
  results$Pv_total <- round(summary(lm_total)$coefficients[2, 4], 4)
  results$RSE_total <- round(sqrt(deviance(lm_total) / df.residual(lm_total)), 4)
  results$tmed_total <- round(mean(ss), 2)
  
  # Análisis por segmentos
  if (length(breakdates) > 0) {
    for (j in 1:length(breakdates)) {
      results[[paste0("year_break_", j)]] <- breakdates[j]
      
      # Definir índices de segmento
      idx_start <- if (j == 1) 1 else bp$breakpoints[j-1] + 1
      idx_end <- bp$breakpoints[j]
      
      # Segmento Pre
      segment_pre <- ss[idx_start:idx_end]
      years_pre <- years_seq[idx_start:idx_end]
      if(length(segment_pre) > 2) { # lm necesita > 2 puntos
        lm_pre <- lm(segment_pre ~ years_pre)
        results[[paste0("P_pre_", j)]] <- round(coef(lm_pre)[2], 4)
        results[[paste0("Pv_pre_", j)]] <- round(summary(lm_pre)$coefficients[2, 4], 4)
        results[[paste0("RSE_pre_", j)]] <- round(sqrt(deviance(lm_pre) / df.residual(lm_pre)), 4)
        results[[paste0("tmed_pre_", j)]] <- round(mean(segment_pre), 2)
      }
      
      # Segmento Post
      idx_post_start <- bp$breakpoints[j]
      idx_post_end <- if (j == length(breakdates)) length(ss) else bp$breakpoints[j+1]
      
      segment_post <- ss[idx_post_start:idx_post_end]
      years_post <- years_seq[idx_post_start:idx_post_end]
      
      if(length(segment_post) > 2) {
        lm_post <- lm(segment_post ~ years_post)
        results[[paste0("P_post_", j)]] <- round(coef(lm_post)[2], 4)
        results[[paste0("Pv_post_", j)]] <- round(summary(lm_post)$coefficients[2, 4], 4)
        results[[paste0("RSE_post_", j)]] <- round(sqrt(deviance(lm_post) / df.residual(lm_post)), 4)
        results[[paste0("tmed_post_", j)]] <- round(mean(segment_post), 2)
      }
    }
  }
  
  return(as.data.frame(results))
}

# --- 2. Cargar y Preprocesar (Corregir años) ---
a <- rast("C:/A_TRABAJO/A_JORGE/ERA5/ERA5_DATA/era5_1940_2025.nc")


# 1. Crear la secuencia de fechas correcta
# Asegúrate de que el número de fechas sea igual a nlyr(a)
n_capas <- nlyr(a)
fechas <- seq(as.Date("1940-01-01"), by = "month", length.out = n_capas)

# 2. Asignar el tiempo usando terra::time
terra::time(a) <- fechas

# 3. Ahora sí funcionará tapp
b <- terra::tapp(a, "years", max)
b <- b - 273.15
data_rast <- terra::rotate(b)
crs(data_rast) <- "epsg:4326"



df_data <- as.data.frame(data_rast, xy = TRUE)
# Limpieza de datos (depende de tu dataset)
df_data <- na.omit(df_data) 

# --- 3. Ejecución Paralela ---
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

tic()
res <- foreach(p = 1:nrow(df_data), .combine = rbind, .packages = c("strucchange", "stats")) %dopar% {
  # Definir años reales de tus datos aquí
  bp_analysis(p, df_data, start_year = 1940, end_year = 2026) 
}
toc()
stopCluster(cl)

final <- cbind(df_data[,c(1,2)], res)
write.csv2(final, "C:/A_TRABAJO/ERA5/RES_ERA5_1940_2025_corregido.csv")