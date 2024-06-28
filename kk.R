library(raster)
library(terra)
library(strucchange)
library(doParallel)
library(foreach)
library(tictoc)

# Cargar los datos
a <- rast("C:/A_TRABAJO/ERA5/ERA5_DATA/era5_1940_2023.nc")
n_layers <- nlyr(a)
n_layers <- seq(1, n_layers, by = 2)
a <- a[[n_layers]]

# Preprocesamiento de datos
b <- terra::tapp(a, "years", max)
b <- b - 273.15
data <- rotate(b)
crs(data) <- "epsg:4326"
data <- as.data.frame(data, xy = TRUE)
data <- data[, -87]


##################
res <- 
  data.frame(
    year_break_1 = NA, year_break_2 = NA, year_break_3 = NA, year_break_4 = NA, year_break_5 = NA,
    n_years = NA, F_sup = NA, p_value = NA,
    P_total = NA, Pv_total = NA, AIC_total = NA, BIC_total = NA, RSE_total = NA,
    P_pre_1 = NA, Pv_pre_1 = NA, AIC_pre_1 = NA, BIC_pre_1 = NA, RSE_pre_1 = NA,
    P_post_1 = NA, Pv_post_1 = NA, AIC_post_1 = NA, BIC_post_1 = NA, RSE_post_1 = NA,
    P_pre_2 = NA, Pv_pre_2 = NA, AIC_pre_2 = NA, BIC_pre_2 = NA, RSE_pre_2 = NA,
    P_post_2 = NA, Pv_post_2 = NA, AIC_post_2 = NA, BIC_post_2 = NA, RSE_post_2 = NA,
    P_pre_3 = NA, Pv_pre_3 = NA, AIC_pre_3 = NA, BIC_pre_3 = NA, RSE_pre_3 = NA,
    P_post_3 = NA, Pv_post_3 = NA, AIC_post_3 = NA, BIC_post_3 = NA, RSE_post_3 = NA,
    P_pre_4 = NA, Pv_pre_4 = NA, AIC_pre_4 = NA, BIC_pre_4 = NA, RSE_pre_4 = NA,
    P_post_4 = NA, Pv_post_4 = NA, AIC_post_4 = NA, BIC_post_4 = NA, RSE_post_4 = NA,
    P_pre_5 = NA, Pv_pre_5 = NA, AIC_pre_5 = NA, BIC_pre_5 = NA, RSE_pre_5 = NA,
    P_post_5 = NA, Pv_post_5 = NA, AIC_post_5 = NA, BIC_post_5 = NA, RSE_post_5 = NA
  )
res <- res[-1,]



# Detectar el número de núcleos disponibles en tu máquina
numCores <- detectCores()

# Registrar el backend de paralelización
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Paralelizar el bucle for
res <- foreach(p = 1:5, .combine = rbind, .packages = c("strucchange", "stats")) %dopar% {
  bp_analysis(p, data)
}

# Detener el clúster
stopCluster(cl)

final <- cbind(data[1:5,c(1,2)], res)


#######################################################
#####FUNCION###################
bp_analysis <- function(i, data) {
  ss <- as.numeric(data[i, -c(1, 2)])
  ss <- ts(ss, start = 1940, end = 2023, frequency = 1)
  
  bp <- strucchange::breakpoints(ss ~ 1, data = ss)
  year_break <- strucchange::breakdates(bp)
  
  # Inicializar data frame de resultados
  results <- list(
    year_break_1 = NA, year_break_2 = NA, year_break_3 = NA, year_break_4 = NA, year_break_5 = NA,
    n_years = NA, F_sup = NA, p_value = NA,
    P_total = NA, Pv_total = NA, AIC_total = NA, BIC_total = NA, RSE_total = NA,
    P_pre_1 = NA, Pv_pre_1 = NA, AIC_pre_1 = NA, BIC_pre_1 = NA, RSE_pre_1 = NA,
    P_post_1 = NA, Pv_post_1 = NA, AIC_post_1 = NA, BIC_post_1 = NA, RSE_post_1 = NA,
    P_pre_2 = NA, Pv_pre_2 = NA, AIC_pre_2 = NA, BIC_pre_2 = NA, RSE_pre_2 = NA,
    P_post_2 = NA, Pv_post_2 = NA, AIC_post_2 = NA, BIC_post_2 = NA, RSE_post_2 = NA,
    P_pre_3 = NA, Pv_pre_3 = NA, AIC_pre_3 = NA, BIC_pre_3 = NA, RSE_pre_3 = NA,
    P_post_3 = NA, Pv_post_3 = NA, AIC_post_3 = NA, BIC_post_3 = NA, RSE_post_3 = NA,
    P_pre_4 = NA, Pv_pre_4 = NA, AIC_pre_4 = NA, BIC_pre_4 = NA, RSE_pre_4 = NA,
    P_post_4 = NA, Pv_post_4 = NA, AIC_post_4 = NA, BIC_post_4 = NA, RSE_post_4 = NA,
    P_pre_5 = NA, Pv_pre_5 = NA, AIC_pre_5 = NA, BIC_pre_5 = NA, RSE_pre_5 = NA,
    P_post_5 = NA, Pv_post_5 = NA, AIC_post_5 = NA, BIC_post_5 = NA, RSE_post_5 = NA
  )
  
  qlr <- strucchange::Fstats(ss ~ 1, data = ss, from = 0.03)
  test <- strucchange::sctest(qlr, type = "supF")
  results$F_sup <- test$statistic
  results$p_value <- test$p.value
  
  # Análisis total
  year_total <- seq(1940, 2023, 1)
  lm_total <- lm(ss ~ year_total)
  results$n_years <- length(year_break)
  results$P_total <- round(lm_total$coefficients[2], 4)
  results$Pv_total <- round(summary(lm_total)$coefficients[2, 4], 4)
  results$AIC_total <- AIC(lm_total)
  results$BIC_total <- BIC(lm_total)
  results$RSE_total <- sqrt(deviance(lm_total) / df.residual(lm_total))
  
  # Análisis por segmentos
  for (j in 1:length(year_break)) {
    if (is.na(year_break[j])) {
      results[[paste0("year_break_", j)]] <- 999
    } else {
      results[[paste0("year_break_", j)]] <- year_break[j]
      pre <- if (j == 1) ss[1:bp$breakpoints[j]] else ss[(bp$breakpoints[j-1] + 1):bp$breakpoints[j]]
      year_pre <- if (j == 1) seq(1940, year_break[j], 1) else seq(year_break[j-1] + 1, year_break[j], 1)
      
      lm_pre <- lm(pre ~ year_pre)
      results[[paste0("P_pre_", j)]] <- round(lm_pre$coefficients[2], 4)
      results[[paste0("Pv_pre_", j)]] <- round(summary(lm_pre)$coefficients[2, 4], 4)
      results[[paste0("AIC_pre_", j)]] <- AIC(lm_pre)
      results[[paste0("BIC_pre_", j)]] <- BIC(lm_pre)
      results[[paste0("RSE_pre_", j)]] <- sqrt(deviance(lm_pre) / df.residual(lm_pre))
      
      post <- if (j == length(year_break)) ss[(bp$breakpoints[j] + 1):length(ss)] else ss[(bp$breakpoints[j] + 1):bp$breakpoints[j+1]]
      year_post <- if (j == length(year_break)) seq(year_break[j] + 1, 2023, 1) else seq(year_break[j] + 1, year_break[j+1], 1)
      
      lm_post <- lm(post ~ year_post)
      results[[paste0("P_post_", j)]] <- round(lm_post$coefficients[2], 4)
      results[[paste0("Pv_post_", j)]] <- round(summary(lm_post)$coefficients[2, 4], 4)
      results[[paste0("AIC_post_", j)]] <- AIC(lm_post)
      results[[paste0("BIC_post_", j)]] <- BIC(lm_post)
      results[[paste0("RSE_post_", j)]] <- sqrt(deviance(lm_post) / df.residual(lm_post))
    }
  }
  
  return(as.data.frame(results))
}

