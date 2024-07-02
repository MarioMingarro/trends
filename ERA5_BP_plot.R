
kk <- final[385341,]
writexl::write_xlsx(kk, "C:/A_TRABAJO/ERA5/kk_data.xlsx")
kk <- data[385341,]
#5 385341
#4
#3
#2 324911
#1 13249
ss <- as.numeric(kk[, -c(1, 2)])
ss <- ts(ss, start = 1940, end = 2023, frequency = 1)


n <- 2023 - 1940 + 1

# Crear la serie de años
years <- 1940:2023

# Convertir la serie temporal en un data frame con años
df <- data.frame(
  year = years,
  value = ss
)

# Calcular los puntos de ruptura
bp <- strucchange::breakpoints(ss ~ 1, data = ss)

# Extraer los puntos de ruptura
breakpoints <- c(1, bp$breakpoints, n)

# Crear un data frame para almacenar los valores ajustados
df_fitted <- data.frame(year = numeric(), fitted_value = numeric())

# Ajustar un modelo lineal en cada segmento
for (i in 1:(length(breakpoints) - 1)) {
  segment <- df %>%
    filter(year >= df$year[breakpoints[i]] & year <= df$year[breakpoints[i + 1]])
  
  if(nrow(segment) > 0) {
    model <- lm(value ~ year, data = segment)
    segment$fitted_value <- predict(model, newdata = segment)
    
    df_fitted <- rbind(df_fitted, segment)
  }
}

# Graficar la serie temporal y las líneas de ajuste con ggplot2
ggplot(df, aes(x = year, y = value)) +
  geom_line(color = "blue") +  # Serie temporal original
  geom_line(data = df_fitted, aes(y = fitted_value), color = "red", size = 1.5) +  # Líneas de ajuste
  theme_minimal() +
  labs(title = "Structural change",
       x = "Year",
       y = "T(ºC)")
