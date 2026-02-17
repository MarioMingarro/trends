library(raster)
library(terra)
library(strucchange)
library(doParallel)
library(foreach)
library(tictoc)


bp_analysis_with_test <- function(i, data_matrix, start_year, end_year) {                
  ss_vals <- data_matrix[i, ]
  ss <- ts(ss_vals, start = start_year, end = end_year, frequency = 1)
  years_seq <- seq(start_year, end_year, 1)
  
  # 1. Detección de rupturas (necesario para breakdates)
  bp <- strucchange::breakpoints(ss ~ 1) 
  breakdates <- strucchange::breakdates(bp)
  n_breaks <- length(breakdates)
  
  # 2. prueba de estabilidad
  qlr <- strucchange::Fstats(ss ~ 1, from = 0.03) #Quandt Likelihood Ratio (QLR)
  test <- strucchange::sctest(qlr, type = "supF")
  
  # 3. Análisis total
  lm_total <- lm(ss ~ years_seq)
  
  results <- list(
    n_breaks = n_breaks,
    F_sup = round(test$statistic, 4),
    p_value = round(test$p.value, 4),
    P_total = round(coef(lm_total)[2], 4),
    Pv_total = round(summary(lm_total)$coefficients[2, 4], 4),
    RSE_total = round(sqrt(deviance(lm_total) / df.residual(lm_total)), 4),
    tmed_total = round(mean(ss), 2)
  )
  
  # 4. Análisis por segmentos (igual que antes)
  if (n_breaks > 0) {
    pts <- c(1, bp$breakpoints, length(ss))
    
    for (j in 1:(length(pts) - 1)) {
      idx_start <- pts[j]
      idx_end <- pts[j+1]
      if (j > 1) idx_start <- idx_start + 1
      
      segment <- ss[idx_start:idx_end]
      years_seg <- years_seq[idx_start:idx_end]
      
      if(length(segment) > 2) {
        lm_seg <- lm(segment ~ years_seg)
        results[[paste0("P_seg_", j)]] <- round(coef(lm_seg)[2], 4)
        results[[paste0("Pv_seg_", j)]] <- round(summary(lm_seg)$coefficients[2, 4], 4)
        results[[paste0("year_start_", j)]] <- years_seq[idx_start]
        results[[paste0("year_end_", j)]] <- years_seq[idx_end]
      }
    }
  }
  
  return(as.data.frame(results))
}
a <- rast("C:/A_TRABAJO/A_JORGE/ERA5/ERA5_DATA/era5_1940_2025.nc")
 ext_test <- ext(-10, 0, 35, 45)
 a <- crop(a, ext_test)

n_capas <- nlyr(a)
fechas <- seq(as.Date("1940-01-01"), by = "month", length.out = n_capas)

terra::time(a) <- fechas


b <- terra::tapp(a, "years", mean)
b <- b - 273.15
data_rast <- terra::rotate(b)
crs(data_rast) <- "epsg:4326"


df_data <- as.data.frame(data_rast, xy = TRUE)
df_data <- na.omit(df_data) 

numCores <- 10
cl <- makeCluster(numCores)
registerDoParallel(cl)

tic()
res <- foreach(p = 1:nrow(df_data), .combine = rbind, .packages = c("strucchange", "stats")) %dopar% {
  bp_analysis(p, df_data, start_year = 1940, end_year = 2025) 
}
toc()
stopCluster(cl)

final <- cbind(df_data[,c(1,2)], res)
write.csv2(final, "C:/A_TRABAJO/ERA5/RES_ERA5_1940_2025_corregido.csv")