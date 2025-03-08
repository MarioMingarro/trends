library(tidyverse)
library(sf)
library(biscale)
library(cowplot)


world_map <- geodata::world(path = tempdir(), resolution = 2)
world_map <- sf::st_as_sf(world_map)
#world_map <- st_transform(world_map, crs = "+proj=robin")
sf_ocean <- rnaturalearth::ne_download(scale = 10, category = "physical", 
                                       type = "ocean",
                                       returnclass = "sf")
#sf_ocean <-st_transform(sf_ocean, crs = "+proj=robin")

res <- read.csv2("C:/A_TRABAJO/A_JORGE/ERA5/ERA5_RESULTS/RESULTADOS_ERA5_1940_2023.csv")
res <- res[,-1]
res[res == 999] <- NA
res <- sf::st_as_sf(res, coords = c("x", "y"), crs = 4326)
res <- st_transform(res, crs = "+proj=robin")

res <- res %>%
  mutate(x = st_coordinates(.)[, "X"], 
         y = st_coordinates(.)[, "Y"],
    Pre_1_C = ifelse(Pv_pre_1 < 0.01, ifelse(P_pre_1 > 0, 1, -1), 0),
    Post_1_C = ifelse(Pv_post_1 < 0.01, ifelse(P_post_1 > 0, 1, -1), 0),
    Post_2_C = ifelse(Pv_post_2 < 0.01, ifelse(P_post_2 > 0, 1, -1), 0),
    Post_3_C = ifelse(Pv_post_3 < 0.01, ifelse(P_post_3 > 0, 1, -1), 0),
    Post_4_C = ifelse(Pv_post_4 < 0.01, ifelse(P_post_4 > 0, 1, -1), 0),
    Post_5_C = ifelse(Pv_post_5 < 0.01, ifelse(P_post_5 > 0, 1, -1), 0))




res_bi <- bi_class(res, x = Pre_1_C, y = Post_1_C, style = "quantile")
#############################################################################################################


res_bi <- bi_class(res, x = Pre_1_C, y = Post_1_C, style = "quantile")

res_bi<- st_transform(res_bi, crs = "+proj=robin")

# Reclasificar bi_class a valores numéricos
res_bi$bi_class_num <- as.numeric(factor(res_bi$bi_class, levels = c(
  "1-1", "1-2", "1-3",
  "2-1", "2-2", "2-3",
  "3-1", "3-2", "3-3"
)))
# Convertir sf a SpatVector de terra
res_bi_terra <- vect(res_bi)

# Crear un SpatRaster vacío con la extensión de tus datos
raster_template <- rast(ext(res_bi_terra), resolution = 0.25) # Ajusta la resolución según sea necesario
crs(raster_template) <- "EPSG:4326"

# Rasterizar los puntos usando bi_class_num como valor
raster_bi <- rasterize(res_bi_terra, raster_template, field = "bi_class_num", fun = "last")
plot(raster_bi)

raster_bi <- project(raster_bi, "+proj=robin")
color_palette <- c(
  "#73ae80",  # 1-1
  "#5a9178",  # 1-2
  "#2a5a5b",  # 1-3
  "#b8d6be",  # 2-1
  "#90b2b3",  # 2-2
  "#567994",  # 2-3
  "#e8e8e8",  # 3-1
  "#b5c0da",  # 3-2
  "#6c83b5"   # 3-3
)
raster_df <- as.data.frame(raster_bi, xy = TRUE)
colnames(raster_df)[3] <- "bi_class_num"

# Plotear el raster con ggplot2
ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = factor(bi_class_num))) +
  geom_sf(data = world_map, fill = "transparent", color = "black",lwd = 0.2) +
  scale_fill_manual(values = color_palette, na.value = "transparent") +
  coord_sf(crs = "+proj=robin") +
  labs(fill = "bi_class") # Añadir leyenda














unique(res_bi$bi_class)
st_write(res_bi,"C:/A_TRABAJO/A_JORGE/ERA5/ERA5_RESULTS/res_bi.shp" )

color_palette <- c(
  "1-1" = "#73ae80",  # Verde claro
  "1-2" = "#5a9178",  # Verde medio
  "1-3" = "#2a5a5b",  # Verde oscuro
  "2-1" = "#b8d6be",  # Verde muy claro
  "2-2" = "#90b2b3",  # Gris claro
  "2-3" = "#567994",  # Azul grisáceo
  "3-1" = "#e8e8e8",  # Gris muy claro
  "3-2" = "#b5c0da",  # Azul claro
  "3-3" = "#6c83b5"   # Azul oscuro
)
ggplot() +
  geom_sf(data = world_map, fill = "gray30", color = "black",lwd = 0.2) +
  geom_sf(data = sf_ocean, fill = "gray80", color = "gray30",lwd = 0.2) +
  geom_sf(data = res_bi, aes(fill = bi_class))+
  scale_fill_manual(values = color_palette, na.value = "transparent") +
  coord_sf(crs = "+proj=robin")


bi_map <- 
  ggplot() +
  geom_sf(data = world_map, fill = "gray30", color = "black",lwd = 0.2) +
  geom_sf(data = sf_ocean, fill = "gray80", color = "gray30",lwd = 0.2) +
  geom_tile(data = res_bi, aes(x = x, y = y, fill = bi_class))+
  bi_scale_fill(pal = "GrPink", dim = 3, guide = F) 

legend <- 
  bi_legend(pal = "DkViolet",
            dim = 3,
            xlab = "Mean ",
            ylab = "Variance ") +
  bi_theme(base_family = "Changa One") +
  theme(rect = element_rect(fill = "grey10"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(size = 21,
                                    color = "#BC7C8F"),
        axis.title.y = element_text(size = 21,
                                    color = "#89A1C8"))

map_legend <- 
  ggdraw() +
  draw_plot(bi_map, 0, 0, 1, 1) +
  draw_plot(legend, 0.15, 0.25, 0.2, 0.2)

ggsave("C:/A_TRABAJO/A_JORGE/ERA5/ERA5_RESULTS/kk.pdf",
       plot = map_legend, width = 15, height = 16.5, device = cairo_pdf)






res_filtrado <- res_filtrado %>%
  mutate(bi_class = bi_class(res_filtrado, x = Pre_1_C, y = Post_1_C, style = "quantile", dim = 3))


legend <- bi_legend(pal = "PinkGrn",
                    xlab = "-1       0       1",
                    ylab = "-1       0       1",
                    size = 10)


ggplot()+
  geom_sf(data = res_filtrado, aes(color = res_filtrado$bi_class))

unique(res_filtrado$bi_class)
ggplot() +
  geom_sf(data = world_map, fill = "gray30", color = "black") +
  geom_sf(data = sf_ocean, fill = "gray80", color = "gray30") +
  geom_point(data = res_filtrado, aes(color = bi_class),size = 0.5) +
  bi_scale_fill(pal = "BlueYl", dim = 3) +
  theme_minimal() +
  annotation_custom(grob = bi_legend, xmin = -18000000, ymin = -8000000)


legend <- bi_legend(pal = "PinkGrn",
                    xlab = "-1       0       1",
                    ylab = "-1       0       1",
                    size = 10)
 
library(biscale)

res <- res %>%
  mutate(bi_class = bi_class(res, x = Pre_1_C, y = Post_1_C, style = "quantile", dim = 3)) # Ajusta dim según necesites
bi_legend <- bi_legend(pal = "BlueYl",
                       dim = 3,
                       xlab = "Pre_1_C",
                       ylab = "Post_1_C",
                       size = 8)
ggplot() +
  geom_sf(data = world_map, fill = "gray30", color = "black") +
  geom_sf(data = sf_ocean, fill = "gray80", color = "gray30") +
  geom_sf(data = res, aes(fill = bi_class), color = NA, size = 0.5) +
  bi_scale_fill(pal = "BlueYl", dim = 3) +
  theme_minimal() +
  annotation_custom(grob = bi_legend, xmin = -18000000, ymin = -8000000) 







theme_set(theme_minimal())




ggplot() +
  geom_sf(data = world_map,
          fill = "gray30",
          color = "black") +
  geom_sf(data = sf_ocean,
          fill = "gray80",
          color = "gray30")+
  geom_point(data = res, aes(x = x, y = y, color = n_years), size = 0.5) +
  theme_minimal()



  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(color = "white", 
                                        fill = "white"),
        plot.background = element_rect(color = "white", 
                                       fill = "white"))
  
  
  geom_tile(data = res, aes(x = x, 
                               y = y, 
                               fill = n_years))

res %>% 
  group_by(n_years) %>% 
  summarize(mean(area....))
world_map = map_data("world") %>% 
  filter(! long > 180)

countries = world_map %>% 
  distinct(region) %>% 
  rowid_to_column()


  ggplot() +
  geom_sf(data = world_map) +
  coord_map("moll")
  