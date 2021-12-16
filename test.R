library(readxl)
library(tidyverse)
library(sp)
library(raster)
library(strucchange)
library(tictoc)
tic()
TXMC <- raster::stack()
for (i in 1980:1990){
  raster <- calc(raster::stack(list.files("B:/DATA/CHELSA/WORLD/TMAX", pattern = paste0(i), full.names = TRUE)), max) # MIN
  TXMC <- raster::stack(TXMC, raster)
}
toc()
names(TXMC) <- paste0("Y_", seq(1979, 2019, by = 1))

mask <- shapefile("C:/GITHUB_REP/butterfly_climate_analysis/Data/Peninsula_Iberica_mask.shp")
random_points <- spsample(mask, n=10, type='random')
data <- raster::extract(TXMC,
                        random_points,
                        df = TRUE)
resultados <- data.frame(year_break = "", 
                         P_pre = "",
                         P_post = "", 
                         P_total = "", 
                         F_sup = "",
                         p.value = "")

for (i in 1:nrow(data)){
  ss <- as.vector(data[i,-1])
  ss <- ts(t(ss),
           start = 1979, 
           end = 2019, 
           frequency = 1)
  #year
  qlr <- Fstats(ss ~ 1, data = ss) #Quandt Likelihood Ratio (QLR)
  bp <- breakpoints(qlr)
  year_break <- breakdates(bp)
  resultados[i,1] <- year_break
  
  
  #trends
  ##pre
  pre <- ss[1:bp$breakpoints]
  year_pre <- seq(1979, year_break, 1)
  lm_pre <- lm(pre ~ year_pre)
  resultados[i,2] <- round(lm_pre$coefficients[2],4)
  
  ##post
  post <- ss[bp$breakpoints:length(ss)]
  year_post <- seq(year_break,2019,1)
  lm_post <- lm(post ~ year_post)
  resultados[i,3] <- round(lm_post$coefficients[2],4)
  
  ##general
  year_total <- seq(1979,2019,1)
  lm_total <- lm(ss ~ year_total)
  resultados[i,4] <- round(lm_total$coefficients[2],4)
  
  # Test the null hypothesis that the annual temperature remains constant over the years
  test <- sctest(qlr, type = "supF")
  resultados[i,5] <- test[1]
  resultados[i,6] <- test[2]
  sa_cusum <- efp(ss ~ 1, data = ss, type = "OLS-CUSUM")
}
kk <- confint(bp, breaks = 1)

###--------------------
sa_cusum <- efp(ss ~ 1, data = ss, type = "OLS-CUSUM")
plot(sa_cusum)
sa_cusum$sigma
my.seg <- segmented(my.lm, seg.Z = ~ year)
summary(sa_cusum)
sa_cusum$coefficients


summary(lm_pre) 

library(forecast)
fit <- tslm(ss[1:24,] ~ trend)
fit$coefficients
#---
decomposedRes <- decompose(ss)
plot(decomposedRes)

#---------------------------------------------------------------#
Datos <- read_excel("C:/GITHUB_REP/trends/Data/Data_test.xlsx", sheet = "Hoja1")
colnames(Datos)

segmented_result <- data.frame(
  Pg ="NULL",
  P1 ="NULL",
  P2 ="NULL",
  Y  ="NULL"
  )

for (i in 1:length(Datos)){
  kk <-  Datos[i,2:117]
  kk <- kk/10
  kk <- as.data.frame(t(kk))
  kk <- mutate(kk, YEAR = seq(1901, 2016, 1))
  colnames(kk) <- c("TMED", "YEAR")
  
  ### SEGMENTED
  my.lm <- lm(TMED ~ YEAR, data = kk)
  my.coef <- coef(my.lm)
  my.seg <- segmented(my.lm, seg.Z = ~ YEAR, psi = list(YEAR = 1950))
  my.fitted <- fitted(my.seg)
  my.model <- data.frame(YEAR = kk$YEAR, TMED = my.fitted)
  
  my.lines <- my.seg$psi[, 2]
  pend <- slope(my.seg)
  Pg <- round(my.lm$coefficients[2],4)
  P1 <- round(pend$YEAR[1,1],4)
  P2 <- round(pend$YEAR[2,1],4)
  Y <- round(my.lines)
  res <- data.frame(Pg,
                    P1,
                    P2,
                    Y )
  segmented_result <- rbind(segmented_result, res)
  
}

### SEGMENTED
my.lm <- lm(TMED ~ YEAR, data = kk)
my.coef <- coef(my.lm)
my.seg <- segmented(my.lm, seg.Z = ~ YEAR, psi = list(YEAR = 1950))
my.fitted <- fitted(my.seg)
my.model <- data.frame(YEAR = kk$YEAR, TMED = my.fitted)

my.lines <- my.seg$psi[, 2]
pend <- slope(my.seg)
Pg <- round(my.lm$coefficients[2],4)
P1 <- round(pend$YEAR[1,1],4)
P2 <- round(pend$YEAR[2,1],4)
Y <- round(my.lines)

ggplot(kk, aes(x = YEAR, y =TMED)) +
  geom_line(col="grey30")+
  geom_point(col="white", alpha=0.8)+
  labs(x = "Año", y= "ºC")+
  geom_smooth(se = FALSE, method = lm)+
  geom_line(data = my.model, aes(x = YEAR, y = TMED), colour = "red", size = 1) +
  geom_vline(xintercept = my.lines, linetype = "dashed", col="black", size=0.5)

plot()
#----------------------

#### strucchange
library(fpp)
library(strucchange)
ss <- ts(kk[,1],frequency=1,start=c(1901,1))
plot(ss)

pp <- breakpoints(ss ~ 1)
summary(pp)
plot(pp)
lines(breakpoints(pp))

ci.nile <- confint(pp)
ci.nile
lines(ci.nile)


sa_bp <- breakpoints(TMED ~ YEAR, data = kk, breaks = 2)  
summary(sa_bp)
kk[90,2]
plot(sa_bp, breaks = 2)
lines(breakpoints(sa_bp))
sa_bp$RSS.table
plot(kk)
lines(breakpoints(sa_bp))
library(mcp)  
model <- list(TMED ~ YEAR)
fit = mcp(model, data = kk)
