library(raster)
library(sf)
library(tidyverse)
library(terra)
library(lubridate) 
library(readxl)
library(strucchange)
library(tictoc)
library(doParallel)
library(foreach)
library(VoCC)
library(tidyterra)

rm(list=(ls()[ls()!="data"]))

# DATA ----
a <- rast("T:/ERA5_DATA/era5_1940_2023.nc")

n_layers <- nlyr(a)

n_layers <- seq(1, n_layers, by = 2)

a <- a[[n_layers]]
rm(n_layers)

tic()
b <- terra::tapp(a, "years", max)
b <- b-273.15
toc() # 15 sec



data <- rotate(b)

crs(data)  <- "epsg:4326"


data <- as.data.frame(data, xy = TRUE)
data <- data[, -87]

data2 <- data
data <- data[490481,]

# BP ANALYSIS ----
tic()
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(i = c(1:nrow(data)), # 1:nrow(data) 1:2500  # 2500:5000 # 5000:7500 # 7500:nrow(data)
               .combine = 'rbind'
) %dopar% {
  ss <- as.vector(data[i,-c(1,2)])
  ss <- ts(t(ss),
           start = 1940,
           end = 2023,
           frequency = 1)
  #year
  qlr <- strucchange::Fstats(ss ~ 1, data = ss, from = 0.05) #Quandt Likelihood Ratio (QLR)
  bp <- strucchange::breakpoints(qlr)
  year_break <- strucchange::breakdates(bp)
  RSS <- bp$RSS
  
  pre <- ss[1:bp$breakpoints]
  year_pre <- seq(1940, year_break, 1)
  lm_pre <- lm(pre ~ year_pre)
  P_pre <- round(lm_pre$coefficients[2],4)
  Pv_pre <- round(summary(lm_pre)$coefficients[2, 4],4)
  AIC_pre <- AIC(lm_pre)
  BIC_pre <- BIC(lm_pre)
  RSE_pre <- sqrt(deviance(lm_pre)/df.residual(lm_pre))
  
  ##post
  post <- ss[bp$breakpoints:length(ss)]
  year_post <- seq(year_break,2023,1)
  lm_post <- lm(post ~ year_post)
  P_post <- round(lm_post$coefficients[2],4)
  Pv_post <- round(summary(lm_post)$coefficients[2, 4],4)
  AIC_post <- AIC(lm_post)
  BIC_post <- BIC(lm_post)
  RSE_post <- sqrt(deviance(lm_post)/df.residual(lm_post))
  
  ##general
  year_total <- seq(1940,2023,1)
  lm_total <- lm(ss ~ year_total)
  P_total <- round(lm_total$coefficients[2],4)
  Pv_total <- round(summary(lm_total)$coefficients[2, 4],4)
  AIC_total <- AIC(lm_total)
  BIC_total <- BIC(lm_total)
  RSE_total <- sqrt(deviance(lm_total)/df.residual(lm_total))
  
  # Test the null hypothesis that the annual temperature remains constant over the years
  test <- strucchange::sctest(qlr, type = "supF")
  F.sup <- test[1]
  p.value <- test[2]
  #sa.cusum <- strucchange::efp(ss ~ 1, data = ss, type = "OLS-CUSUM")
  data.frame(year_break, RSS, F.sup, p.value, 
             P_total, Pv_total, AIC_total, BIC_total, RSE_total,
             P_pre, Pv_pre, AIC_pre, BIC_pre,RSE_pre,
             P_post, Pv_post,AIC_post, BIC_post, RSE_post)
}

parallel::stopCluster(cl = my.cluster)
toc()

res <- cbind(data, res)
resultados <- res


write.csv2(res,"D:/ERA_5_2024/resultados_ERA5_max2.csv" )


sig <- dplyr::filter(resultados, p.value <= 0.01)
sig <- dplyr::filter(resultados, resultados$ <= 0.01)

# PLOTTING ----

#year_break
#RSS
#statistic
#p.value
#P_total
#AIC_total
#BIC_total
#RSE_total
#P_pre
#AIC_pre
#BIC_pre
#RSE_pre
#P_post
#AIC_post
#BIC_post
#RSE_post



## Raster creation----
s <- c("year_break", "statistic", "p.value", "P_total", "Pv_total", "RSE_total", "P_pre", "RSE_pre", "Pv_pre","P_post", "RSE_post", "Pv_post")
f <- c("x", "y", "year_break", "statistic", "p.value", "P_total", "Pv_total", "RSE_total", "P_pre", "RSE_pre", "Pv_pre","P_post", "RSE_post", "Pv_post")

resultados_land <-res
resultados_land <- subset(sig, select=f)


for (i in 1:length(s)){
  ## Land
  r <- data.frame(x=resultados_land[,1], y=resultados_land[,2], resultados_land[i+2])
  r <- rast(r, type="xyz")
  crs(r)  <- "epsg:4326"
  
 writeRaster(r,paste0("D:/ERA_5_2024/MAPS/SIG/", s[i] ,"_max.tif" ))

r <- project(r,"+proj=hatano", mask = TRUE)


# Discretize for better plotting after projection
g <- st_graticule(ndiscr = 1000) %>%  st_transform(st_crs(r))

border <- st_graticule() %>% 
  st_bbox()%>%
  st_as_sfc()%>%
  st_transform(3857)%>%
  st_segmentize(500000)%>%
  st_transform(st_crs(r))%>%
  st_cast("POLYGON")

# Labels
labels_x_init <- g %>%
  filter(type == "N") %>%
  mutate(lab = paste0(degree, "°"))

labels_x <- st_as_sf(st_drop_geometry(labels_x_init), lwgeom::st_startpoint(labels_x_init))


labels_y_init <- g %>%
  filter(type == "E") %>%
  mutate(lab = paste0(degree, "°"))

labels_y <- st_as_sf(st_drop_geometry(labels_y_init), lwgeom::st_startpoint(labels_y_init))



# Plot
print(ggplot() +
  geom_sf(data = border, fill = "azure", color = "lightgray", linewidth = .5) +
  geom_sf(data = g, color = "lightgray") +
  geom_spatraster(data = r) +
  scale_fill_whitebox_c(palette = "viridi") +
  geom_sf_text(data = labels_x, aes(label = lab), nudge_x = -1000000, size = 3) +
  geom_sf_text(data = labels_y, aes(label = lab), nudge_y = -1000000, size = 3) +
  theme_void() +
  labs(x = "", y = "", fill = paste0(s[i])))

ggsave(paste0("D:/ERA_5_2024/MAPS/SIG/", s[i] ,"_max.jpeg" ), dpi = 600, width = 30, height = 30, units = "cm")

}
