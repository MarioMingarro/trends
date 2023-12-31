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

a <- tapp(a, "years", mean)
a <- a-273.15
b <- tapp(b, "years", mean)
b <- b-273.15

data <- c(a,b)
data <- rotate(data)

crs(data)  <- "epsg:4326"


land <- vect("A:/ERA5_TEST/ne_10m_land.shp")
ocean <- vect("A:/ERA5_TEST/ne_10m_ocean.shp")
#eu <- vect("A:/ERA5_TEST/europe_r.shp")

terra::ext(data) <- terra::ext(ocean)

crs(data)  <- "epsg:4326"

data_L <- terra::mask(data, land)



## Climate velocity ----


OST <- data
crs(OST)  <- "epsg:4326"

vel_result <- as.data.frame(data, xy = TRUE)[,1:2]

for(i in 1:82){ #:nlyr(d)
  
  ia <- i+1
  OST1 <- data[[i]]
  OST2 <- data[[ia]]
  OST <- c(OST1, OST2)
  OST <- raster::stack(OST)
  
  vt <- tempTrend(OST,
                  th = nlayers(OST))
  
  vg <- spatGrad(OST,
                 th = 0.1,
                 projected = FALSE)
  
  gv <- gVoCC(vt, 
              vg)
  
  vel <- gv[[1]]
  kk2 <- as.data.frame(vel, xy = T)
  vel_result <- cbind(vel_result, kk2[,3])
  
}

colnames(vel_result) <-  c("x", "y", names(data)[2:83])

write.csv2(acc_result,"A:/ERA5_RESULTS/acc_result.csv")


acc_result <- as.data.frame(data, xy = TRUE)[,1:2]
for(i in 3:83){
  ia <- i+1
  acc <- (vel_result[,i]-vel_result[ia])/1
  acc_result <- cbind(acc_result, acc)
}

kk <- rowMeans(vel_result[3:84])
vel_result <- cbind(vel_result, kk)

x2.df <- data.frame(x=vel_result[,1], y=vel_result[,2], vel_result$)
x2 <- rast(x2.df, type="xyz")
crs(x2)  <- "epsg:4326"
plot(x2)

x2 <- terra::mask(x2, ocean)
terra::ext(x2) <- terra::ext(ocean)

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


plot(b)
