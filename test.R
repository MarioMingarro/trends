library(readxl)
library(tidyverse)
library(segmented)

Datos <- read_excel("C:/GITHUB_REP/trends/Data/Data_test.xlsx", sheet = "Hoja1")
colnames(Datos)
kk <-  Datos[1,2:117]
kk <- kk/10
kk <- as.data.frame(t(kk))
kk <- mutate(kk, YEAR = seq(1901, 2016, 1))
colnames(kk) <- c("TMED", "YEAR")


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
  geom_line(data = my.model, aes(x = Distance, y = Elevation), colour = "red", size = 1) +
  geom_vline(xintercept = my.lines, linetype = "dashed", col="black", size=0.5)

plot()
#----------------------


  
  
library(mcp)  
model <- list(TMED ~ YEAR)
fit = mcp(model, data = kk)
