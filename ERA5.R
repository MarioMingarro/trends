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

# DATA ----

a <- rast("B:/A_DATA/ERA_5/ERA_5_1940_1987.nc")
b <- rast("B:/A_DATA/ERA_5/ERA_5_1988_2023.nc")
b <- b[[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,
          77,79,81,83,85,87,89,91,93,95,97,99,101,103,105,107,109,111,113,115,117,119,121,123,125,127,129,131,133,135,
          137,139,141,143,145,147,149,151,153,155,157,159,161,163,165,167,169,171,173,175,177,179,181,183,185,187,189,
          191,193,195,197,199,201,203,205,207,209,211,213,215,217,219,221,223,225,227,229,231,233,235,237,239,241,243,
          245,247,249,251,253,255,257,259,261,263,265,267,269,271,273,275,277,279,281,283,285,287,289,291,293,295,297,
          299,301,303,305,307,309,311,313,315,317,319,321,323,325,327,329,331,333,335,337,339,341,343,345,347,349,351,
          353,355,357,359,361,363,365,367,369,371,373,375,377,379,381,383,385,387,389,391,393,395,397,399,401,403,405,
          407,409,411,413,415,417,419,421,423,425,427,429,431,433,435,437,439,441,443,445,447,449,451,453,455,457,459,
          461,463,465,467,469,471,473,475,477,479,481,483,485,487,489,491,493,495,497,499,501,503,505,507,509,511,513,
          515,517,519,521,523,525,527,529,531,533,535,537,539,541,543,545,547,549,551,553,555,557,559,561,563,565,567,
          569,571,573,575,577,579,581,583,585,587,589,591,593,595,597,599,601,603,605,607,609,611,613,615,617,619,621,
          623,625,627,629,631,633,635,637,639,641,643,645,647,649,651,653,655,657,659,661,663,665,667,669,671,673,675,
          677,679,681,683,685,687,689,691,693,695,697,699,701,703,705,707,709,711,713,715,717,719,721,723,725,727,729,
          731,733,735,737,739,741,743,745,747,749,751,753,755,757,759,761,763,765,767,769,771,773,775,777,779,781,783,
          785,787,789,791,793,795,797,799,801,803,805,807,809,811,813,815,817,819,821,823,825,827,829,831,833,835,837,
          839)]]

x <- tapp(a, "years", mean)
x <- x-273.15
y <- tapp(b, "years", mean)
y <- y-273.15


data <- c(x,y)
crs(data)  <- "epsg:4326"

land <- vect("A:/ERA5_TEST/ne_10m_land.shp")
ocean <- vect("A:/ERA5_TEST/ne_10m_ocean.shp")
#eu <- vect("A:/ERA5_TEST/europe_r.shp")

terra::ext(data) <- terra::ext(ocean)

data <- terra::mask(data, land)


data <- as.data.frame(data, xy = TRUE)



# ANALYSIS ----


tic()
n.cores <- parallel::detectCores() - 2
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
           end = 2022,
           frequency = 1)
  #year
  qlr <- strucchange::Fstats(ss ~ 1, data = ss) #Quandt Likelihood Ratio (QLR)
  bp <- strucchange::breakpoints(qlr)
  year_break <- strucchange::breakdates(bp)
  RSS <- bp$RSS
  
  pre <- ss[1:bp$breakpoints]
  year_pre <- seq(1940, year_break, 1)
  lm_pre <- lm(pre ~ year_pre)
  P_pre <- round(lm_pre$coefficients[2],4)
  AIC_pre <- AIC(lm_pre)
  BIC_pre <- BIC(lm_pre)
  RSE_pre <- sqrt(deviance(lm_pre)/df.residual(lm_pre))
  
  ##post
  post <- ss[bp$breakpoints:length(ss)]
  year_post <- seq(year_break,2022,1)
  lm_post <- lm(post ~ year_post)
  P_post <- round(lm_post$coefficients[2],4)
  AIC_post <- AIC(lm_post)
  BIC_post <- BIC(lm_post)
  RSE_post <- sqrt(deviance(lm_post)/df.residual(lm_post))
  
  ##general
  year_total <- seq(1940,2022,1)
  lm_total <- lm(ss ~ year_total)
  P_total <- round(lm_total$coefficients[2],4)
  AIC_total <- AIC(lm_total)
  BIC_total <- BIC(lm_total)
  RSE_total <- sqrt(deviance(lm_total)/df.residual(lm_total))
  
  # Test the null hypothesis that the annual temperature remains constant over the years
  test <- strucchange::sctest(qlr, type = "supF")
  F.sup <- test[1]
  p.value <- test[2]
  #sa.cusum <- strucchange::efp(ss ~ 1, data = ss, type = "OLS-CUSUM")
  data.frame(year_break, RSS, F.sup, p.value, 
             P_total, AIC_total, BIC_total, RSE_total,
             P_pre, AIC_pre, BIC_pre,RSE_pre,
             P_post,AIC_post, BIC_post, RSE_post)
}

parallel::stopCluster(cl = my.cluster)
toc()
resultados_ocean <- cbind(data, res)

unique(resultados$year_break)
write.csv2(resultados,"A:/ERA5_TEST/ocean_results.csv" )
resultados_land_t <- read.csv2("A:/ERA5_TEST/land_results.csv")
x2.df <- data.frame(x=resultados_ocean[,1], y=resultados_ocean[,2], resultados_ocean$P_pre)
x2 <- rast(x2.df, type="xyz")
crs(x2)  <- "epsg:4326"


# plotting ----
r <- project(x2,"+proj=hatano", mask = TRUE)


# Discretize for better plotting after projection
g <- st_graticule(ndiscr = 500) |> st_transform(st_crs(r))
border <- st_graticule() |>
  st_bbox() |>
  st_as_sfc() |>
  st_transform(3857) |>
  st_segmentize(500000) |>
  st_transform(st_crs(r)) |>
  st_cast("POLYGON")

# Get label placement,
# This is the hardest part
library(dplyr)
labels_x_init <- g %>%
  filter(type == "N") %>%
  mutate(lab = paste0(degree, "°"))

labels_x <- st_as_sf(st_drop_geometry(labels_x_init), lwgeom::st_startpoint(labels_x_init))


labels_y_init <- g %>%
  filter(type == "E") %>%
  mutate(lab = paste0(degree, "°"))

labels_y <- st_as_sf(st_drop_geometry(labels_y_init), lwgeom::st_startpoint(labels_y_init))


# Plot
ggplot() +
  geom_sf(data = border, fill = "azure", color = "lightgray", linewidth = 1) +
  geom_sf(data = g, color = "lightgray") +
  geom_sf(data=land, color = "lightgray")+
  geom_spatraster(data = r) +
  scale_fill_whitebox_c(palette = "viridi") +
  geom_sf_text(data = labels_x, aes(label = lab), nudge_x = -1000000, size = 3) +
  geom_sf_text(data = labels_y, aes(label = lab), nudge_y = -1000000, size = 3) +
  theme_void() +
  labs(x = "", y = "", fill = "Temp")




# Discretize for better plotting after projection
g <- st_graticule(ndiscr = 500) |> st_transform(st_crs(r))
border <- st_graticule() |>
  st_bbox() |>
  st_as_sfc() |>
  st_transform(3857) |>
  st_segmentize(500000) |>
  st_transform(st_crs(r)) |>
  st_cast("POLYGON")

# Get label placement,
# This is the hardest part
library(dplyr)
library(tidyterra)
labels_x_init <- g %>%
  filter(type == "N") %>%
  mutate(lab = paste0(degree, "°"))

labels_x <- st_as_sf(st_drop_geometry(labels_x_init), lwgeom::st_startpoint(labels_x_init))


labels_y_init <- g %>%
  filter(type == "E") %>%
  mutate(lab = paste0(degree, "°"))

labels_y <- st_as_sf(st_drop_geometry(labels_y_init), lwgeom::st_startpoint(labels_y_init))


# Plot
ggplot() +
  geom_sf(data = border, fill = "azure", color = "lightgray", linewidth = 1) +
  geom_sf(data = g, color = "lightgray") +
  tidyterra::geom_spatraster(data = r) +
  scale_fill_whitebox_c(palette = "viridi") +
  geom_sf_text(data = labels_x, aes(label = lab), nudge_x = -1000000, size = 3) +
  geom_sf_text(data = labels_y, aes(label = lab), nudge_y = -1000000, size = 3) +
  theme_void() +
  labs(x = "", y = "", fill = "Temp")

