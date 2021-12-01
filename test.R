library(readxl)
library(tidyverse)
library(segmented)

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
