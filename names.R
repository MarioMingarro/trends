library(dplyr)
library(sf) # Asumiendo que 'geometry' es un objeto sf

# 1. Renombrar todas las columnas homogeneizando pr/pre y pst/post
res <- res %>%
  rename(
    # --- Puntos de ruptura y estadísticos globales ---
    year_break_1 = yr_br_1,
    year_break_2 = yr_br_2,
    year_break_3 = yr_br_3,
    year_break_4 = yr_br_4,
    year_break_5 = yr_br_5,
    p_value_Fsup = p_value,
    
    # --- Tendencia Total ---
    slope_total = P_total,
    p_value_total = Pv_totl,
    rse_total = RSE_ttl,
    
    # --- Ruptura 1 ---
    p_value_pre_1 = Pv_pr_1,
    rse_pre_1 = RSE_pr_1,
    slope_post_1 = P_pst_1,
    p_value_post_1 = Pv_ps_1,
    rse_post_1 = RSE_ps_1,
    
    # --- Ruptura 2 ---
    slope_pre_2 = P_pre_2,
    p_value_pre_2 = Pv_pr_2,
    rse_pre_2 = RSE_pr_2,
    slope_post_2 = P_pst_2,
    p_value_post_2 = Pv_ps_2,
    rse_post_2 = RSE_ps_2,
    
    # --- Ruptura 3 ---
    slope_pre_3 = P_pre_3,
    p_value_pre_3 = Pv_pr_3,
    rse_pre_3 = RSE_pr_3,
    slope_post_3 = P_pst_3,
    p_value_post_3 = Pv_ps_3,
    rse_post_3 = RSE_ps_3,
    
    # --- Ruptura 4 ---
    slope_pre_4 = P_pre_4,
    p_value_pre_4 = Pv_pr_4,
    rse_pre_4 = RSE_pr_4,
    slope_post_4 = P_pst_4,
    p_value_post_4 = Pv_ps_4,
    rse_post_4 = RSE_ps_4,
    
    # --- Ruptura 5 ---
    slope_pre_5 = P_pre_5,
    p_value_pre_5 = Pv_pr_5,
    rse_pre_5 = RSE_pr_5,
    slope_post_5 = P_pst_5,
    p_value_post_5 = Pv_ps_5,
    rse_post_5 = RSE_ps_5,
    
    # --- Temperaturas medias ---
    temp_total = tmd_ttl,
    temp_pre_1 = tmd_pr_1,
    temp_post_1 = tmd_ps_1,
    temp_pre_2 = tmd_pr_2,
    temp_post_2 = tmd_ps_2,
    temp_pre_3 = tmd_pr_3,
    temp_post_3 = tmd_ps_3,
    temp_pre_4 = tmd_pr_4,
    temp_post_4 = tmd_ps_4,
    temp_pre_5 = tmd_pr_5,
    temp_post_5 = tmd_ps_5,
    
    # --- Clasificaciones y Cobertura ---
    area_land = are____,
    slope_bi_pre_1 = Pre_1_C,
    slope_bi_post_1 = Pst_1_C,
    slope_bi_post_2 = Pst_2_C,
    slope_bi_post_3 = Pst_3_C,
    slope_bi_post_4 = Pst_4_C,
    slope_bi_post_5 = Pst_5_C,
    class_change_1 = b_P_1_C_1,
    class_change_2 = b_P_1_C_2,
    class_change_3 = b_P_1_C_3,
    class_change_4 = b_P_1_C_4,
    class_change_5 = b_P_1_C_5
  )

# 2. Extraer coordenadas y añadirlas al principio
# Esto asume que 'res' tiene una columna 'geometry' de tipo sf
res_coords <- res %>%
  mutate(
    lon = st_coordinates(geometry)[, 1],
    lat = st_coordinates(geometry)[, 2]
  ) %>%
  select(lon, lat, everything()) %>% # Poner lon/lat al principio
  select(-geometry) # Quitar la columna geométrica para CSV

# 3. Guardar como archivo CSV
write.csv2(res_coords, "C:/A_TRABAJO/A_JORGE/ERA5/ERA5_RESULTS/resultados_ERA5_1940_2023_nuevo.csv", row.names = TRUE)

