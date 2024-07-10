library(raster)
library(terra)
library(tidyverse)
library(rnaturalearth)
library(sf)
library(ggpubr)
library(gridExtra)

a <- rast("C:/A_TRABAJO/ERA5/ERA5_DATA/era5_1940_2023.nc")
n_layers <- nlyr(a)
n_layers <- seq(1, n_layers, by = 2)
a <- a[[n_layers]]

# Preprocesamiento de datos
b <- terra::tapp(a, "years", max)
b <- b - 273.15
data <- terra::rotate(b)
crs(data) <- "epsg:4326"
data <- as.data.frame(data, xy = TRUE)
data <- data[, -87]
rm(a, b, n_layers)

final <- read.csv2("C:/A_TRABAJO/ERA5/RESULT_ERA5_1940_2023.csv")
final <- final[,-1]


# Seleccionar caso especifico
fila <- 189749
ss <-  data[fila,]
tt <- final[fila,]

#5 385341
#4 365220 829429 586578 566396
#3 380883 750951 622977
#2 32027 189749 225096
#1 13249  19698 198000
#0 396964
###### -----
ss <- as.numeric(ss[, -c(1, 2)])
ss <- ts(ss, start = 1940, end = 2023, frequency = 1)

n <- 2023 - 1940 + 1
# Crear la serie de años
years <- 1940:2023

# Convertir la serie temporal en un data frame con años
df <- data.frame(
  year = years,
  value = ss)

# Calcular los puntos de ruptura
bp <- strucchange::breakpoints(ss ~ 1, data = ss)

# Extraer los puntos de ruptura
breakpoints <- c(1, bp$breakpoints, n)

# Crear un data frame para almacenar los valores ajustados
df_fitted <- data.frame(year = numeric(), fitted_value = numeric(), segment = numeric())

# Ajustar un modelo lineal en cada segmento
for (i in 1:(length(breakpoints) - 1)) {
  segment <- df %>%
    filter(year >= df$year[breakpoints[i]] & year <= df$year[breakpoints[i + 1]])
  
  if (nrow(segment) > 0) {
    model <- lm(value ~ year, data = segment)
    segment$fitted_value <- predict(model, newdata = segment)
    segment$segment <- i  # Añadir una columna para identificar el segmento
    
    df_fitted <- rbind(df_fitted, segment)
  }
}

# Graficar la serie temporal y las líneas de ajuste con ggplot2
trend_plot <- ggplot(df, aes(x = year, y = value)) +
  geom_line(color = "blue") +  # Serie temporal original
  geom_line(data = df_fitted, aes(y = fitted_value, group = segment), color = "red", size = 1.5, alpha = 0.5) +  # Líneas de ajuste por segmento
  theme_minimal() +
  labs(title = "Structural change",
       x = "Year",
       y = "T(ºC)")


# Map
world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")
tt[,1:2] <- as.numeric(tt[,1:2])
# Convertir el data frame del punto a un objeto sf
tt_sf <- st_as_sf(tt, coords = c("x", "y"), crs = 4326)  # Suponiendo que las coordenadas están en WGS84 (EPSG:4326)

# Transformar el sistema de coordenadas del punto al mismo que el mapa del mundo
tt_sf_transformed <- st_transform(tt_sf, crs = st_crs(world))

# Crear el mapa y agregar el punto transformado
map <- ggplot() +
  geom_sf(data = world) +
  theme_minimal() +
  geom_sf(data = tt_sf_transformed, col = "red", fill = "red", shape = 21, size = 2) +
  labs(x = NULL, y = NULL)



#Tabla
tt <- final[fila,]
tt <- tt %>%
  pivot_longer(
    cols = matches("year_break_\\d+|P_pre_\\d+|Pv_pre_\\d+|RSE_pre_\\d+|P_post_\\d+|Pv_post_\\d+|RSE_post_\\d+"),
    names_to = c("metric", "year_break"),
    names_pattern = "(year_break|P_pre|Pv_pre|RSE_pre|P_post|Pv_post|RSE_post)_(\\d+)"
  )

# Ordenar el dataframe por x, y, year_break
tt <- tt %>%
  arrange(x, y, year_break)
tt <- tt[, c(9,11)]
tt <- tt[-c(9,10,11,16,17,18,23,24,25,30,31,32),]

# Identificar los bloques por year_break
tt$group <- cumsum(tt$metric == "year_break")

# Filtrar y pivotear los datos
tt <- tt %>%
  filter(!is.na(value)) %>%
  pivot_wider(names_from = metric, values_from = value, values_fill = NA) %>%
  select(-group)

# Ver el resultado
tt <- as.data.frame(t(tt))
names(tt) <- gsub(x = names(tt), pattern = "V", replacement = "BP_")  

tt[2:7,] <- round(tt[2:7,] ,4)
tk <- ttheme_minimal(
  colhead=list(fg_params=list(col="gray30", fontface=4L)),
  rowhead=list(fg_params=list(col="gray50", fontface=3L)))
table <- tableGrob(tt, theme = tk) 
# BIC

bic <- summary(bp)
bic <- as.data.frame(t(bic$RSS))
bic <- dplyr::mutate(bic, "breaks" = rownames(bic))


# bic <- ggplot(bic, aes(x = breaks)) +
#   geom_line(aes(y = BIC, color = "BIC")) +
#   geom_point(aes(y = BIC, color = "BIC")) +
#   geom_line(aes(y = RSS * 100, color = "RSS")) +  # Escalamos RSS para mejor visualización
#   geom_point(aes(y = RSS * 100, color = "RSS")) +  # Escalamos RSS para mejor visualización
#   scale_y_continuous(sec.axis = sec_axis(~./100, name = "RSS")) +  # Agregar segundo eje y para RSS
#   labs(
#     x = "Número de Puntos de Ruptura",
#     y = "BIC",
#     color = "Métrica"
#   ) +
#   scale_color_manual(values = c("BIC" = "blue", "RSS" = "red")) + 
#   theme_minimal()  # Estilo del gráfico


table_bic <- tableGrob(round(bic[,c(2,1)],4), theme = tk) 

###### -----
# Todos
ggarrange(trend_plot, table_bic, map,
          table)




# Histogram -----
final <- read.csv2("C:/A_TRABAJO/ERA5/RESULT_ERA5_1940_2023.csv")
final <- final[,-1]
sig <- dplyr::filter(final, p_value <= 0.05)
sig$n_years[sig$n_years == "NA"] <- 0
hist(sig$n_years)
View(final)


p1 <- hist(final$n_years)
p2 <- hist(sig$n_years)
plot( p1, col=rgb(0,0,1,1/4))
plot( p2, col=rgb(1,0,0,1/4), add=T)

options(scipen = 999)
library(gridExtra)
ggplot() + 
  geom_histogram(data = final, aes(x = n_years), alpha = 0.3, fill = "red") +
  geom_histogram(data = sig, aes(x = n_years), alpha = 0.3, fill ="blue") +
  labs(x = "Number of breaking points", y= "Number of pixels")+
  theme_bw()+
  annotation_custom(tableGrob(kk), xmin=3, xmax=5, ymin=400000, ymax=600000)

a <- final %>% 
  group_by(n_years) %>% 
  summarise(n = sum(n_years))

b <- sig %>% 
  group_by(n_years) %>% 
  summarise(n = sum(n_years))

kk <- cbind(a,b)
kk <- kk[,-3]
colnames(kk) <- c("n", "all", "sig")
kk <- mutate(kk, dif = kk$all-kk$sig)
kk <- kk[,-1]


library(ggplot2)
library(gridExtra)
set.seed(1)
mydata <- data.frame(a=1:50, b=rnorm(50))
mytable <- cbind(sites=c("site 1","site 2","site 3","site 4"),mydata[10:13,])
k <- ggplot(mydata,aes(x=a,y=b)) + 
  geom_point(colour="blue") + 
  geom_point(data=mydata[10:13, ], aes(x=a, y=b), colour="red", size=5) + 
  annotation_custom(tableGrob(mytable), xmin=35, xmax=50, ymin=-2.5, ymax=-1)


## TESTING -----
fila <- 19698
test <-  data[fila,]
# Eliminar las primeras dos columnas del dataframe 'test'
test <- test[,-c(1,2)]

# Seleccionar las primeras 19 columnas del dataframe 'test'
test2 <- test[,40:52]

# Crear un vector con los años de 1940 a 1958
year <- seq(1979, 1991, 1)

# Asegurarse de que 'test2' tiene los nombres de fila adecuados
rownames(test2) <- NULL

# Unir el dataframe 'test2' con el vector 'year' como una nueva columna
test2 <- as.data.frame(t(rbind(test2, year)))

k <- lm(V1~V2, data = test2)
summary(k)
writexl::write_xlsx(test, "C:/A_TRABAJO/ERA5/kkdvk2.xlsx")

plot(bp)


## plot segements 

# Obtener los puntos de quiebre
breakpoints <- bp$breakpoints

# Crear un dataframe a partir de la serie temporal
ts_data <- data.frame(
  time = time(ss),
  value = as.numeric(ss)
)

# Agregar un indicador de segmento a los datos
ts_data$segment <- cut(1:nrow(ts_data), breaks = c(0, breakpoints, nrow(ts_data)), labels = FALSE, include.lowest = TRUE)

ggplot(ts_data, aes(x = time, y = value)) +
  geom_line(aes(color = factor(segment))) +
  geom_vline(xintercept = ts_data$time[breakpoints], linetype = "dashed") +
  labs(title = "Segmentos de la Serie Temporal con Puntos de Quiebre",
       x = "Tiempo",
       y = "Valor",
       color = "Segmento") +
  theme_minimal()
