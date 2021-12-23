closeAllConnections()
rm(list=(ls()[ls()!="K"]))
gc(reset=TRUE)
source("Dependencies/Functions.R")

# Maximun temperatura of the hottest month
tic()
TXMC <- raster::stack()
for (i in 1901:2016){
  raster <- raster::aggregate(raster::stack(list.files("B:/DATA/CHELSA/SPAIN/TMAX", pattern = paste0(i), full.names = TRUE)), fact=10, fun=mean)
  raster <- reclassify(raster, c(-Inf, -999, NA))
  raster <- calc(raster, max)
  TXMC <- raster::stack(TXMC, raster)
}

names(TXMC) <- paste0("Y_", seq(1901, 2016, by = 1))

long_lat <- rasterToPoints(TXMC[[1]], spatial = TRUE)
data <- raster::extract(TXMC,
                        long_lat,
                        df = TRUE)

write.csv(data, "C:/GITHUB_REP/trends/spain_tmax.csv")

rm(raster)
rm(TXMC)
rm(long_lat)
gc(reset=TRUE)


tic()
data<- read_csv("spain_tmax.csv")
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(i = 7500:nrow(data), # 1:2500  # 2500:5000 # 5000:7500 # 7500:nrow(data)
               .combine = 'rbind'
) %dopar% {
  ss <- as.vector(data[i,-1])
  ss <- ts(t(ss),
           start = 1901,
           end = 2016,
           frequency = 1)
  #year
  qlr <- strucchange::Fstats(ss ~ 1, data = ss) #Quandt Likelihood Ratio (QLR)
  bp <- strucchange::breakpoints(qlr)
  year_break <- strucchange::breakdates(bp)
  
  pre <- ss[1:bp$breakpoints]
  year_pre <- seq(1901, year_break, 1)
  lm_pre <- lm(pre ~ year_pre)
  P_pre <- round(lm_pre$coefficients[2],4)
  
  
  ##post
  post <- ss[bp$breakpoints:length(ss)]
  year_post <- seq(year_break,2016,1)
  lm_post <- lm(post ~ year_post)
  P_post <- round(lm_post$coefficients[2],4)
  
  ##general
  year_total <- seq(1901,2016,1)
  lm_total <- lm(ss ~ year_total)
  P_total <- round(lm_total$coefficients[2],4)
  
  # Test the null hypothesis that the annual temperature remains constant over the years
  test <- strucchange::sctest(qlr, type = "supF")
  F.sup <- test[1]
  p.value <- test[2]
  #sa.cusum <- strucchange::efp(ss ~ 1, data = ss, type = "OLS-CUSUM")
  data.frame(year_break, P_pre, P_post, P_total, F.sup, p.value)
}

parallel::stopCluster(cl = my.cluster)
toc()
#resultados <- res
resultados <- rbind(resultados, res)
resultados <- resultados[1:nrow(data),]


# Plot
long_lat2 <- as.data.frame(rasterToPoints(raster::aggregate(raster("B:/DATA/CHELSA/SPAIN/TMAX/CHELSAcruts_tmax_1_1902_V.1.0.tif"), 10)))
long_lat2 <- long_lat2[1:nrow(resultados),-3]
long_lat2 <- cbind(long_lat2, id= rownames(long_lat2))
resultados_2 <- cbind(resultados, long_lat2$id)
kk <- cbind(resultados_2, long_lat2)
kk <- data.frame(x = kk$x, y = kk$y, z = kk$year_break)
ggplot(kk, aes(x = x, y = y, col=z))+
  geom_point()+
  scale_colour_viridis_c()
