# =====================================================
# 05_multivariate_analysis.R
#
# Multivariate analysis of climate variability mode
# dominance patterns.
#
# Methods:
#   - Principal component analysis (PCA)
#   - K-means clustering
#   - Spatial visualization
#
# Input:
#   dominance_ratio.csv
#
# Output:
#   cluster_profile.csv
#   climate_type_distribution.csv
#   PCA results
#
# =====================================================



# =====================================================
# 1. Load packages
# =====================================================


library(dplyr)
library(FactoMineR)
library(factoextra)

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)



# =====================================================
# 2. Load dominance contribution dataset
# =====================================================


dominance_ratio <- read.csv(
  
  "dominance_ratio.csv",
  
  stringsAsFactors = FALSE
  
)



# =====================================================
# 3. Prepare PCA matrix
# =====================================================


ratio_mat <- dominance_ratio %>%
  
  select(
    ends_with("_ratio")
  ) %>%
  
  as.matrix()



# =====================================================
# 4. Principal component analysis
# =====================================================


pca_res <-
  
  PCA(
    
    ratio_mat,
    
    scale.unit = TRUE,
    
    graph = FALSE
    
  )



# PCA visualization

pca_plot <-
  
  fviz_pca_biplot(
    
    pca_res,
    
    repel = TRUE
    
  )



print(pca_plot)



# Save PCA scores

pca_scores <-
  
  as.data.frame(
    pca_res$ind$coord
  )


write.csv(
  
  pca_scores,
  
  "PCA_scores.csv",
  
  row.names = FALSE
  
)



# =====================================================
# 5. K-means clustering
# =====================================================


set.seed(123)


km_res <-
  
  kmeans(
    
    ratio_mat,
    
    centers = 3,
    
    nstart = 50
    
  )



dominance_ratio$cluster <-
  
  km_res$cluster




# =====================================================
# 6. Characterize climate dominance types
# =====================================================



cluster_profile <-
  
  dominance_ratio %>%
  
  group_by(cluster) %>%
  
  summarise(
    
    ENSO =
      mean(
        ENSO_intensity_ratio,
        na.rm = TRUE
      ),
    
    
    PDO =
      mean(
        PDO_intensity_ratio,
        na.rm = TRUE
      ),
    
    
    NAO =
      mean(
        NAO_intensity_ratio,
        na.rm = TRUE
      ),
    
    
    AO =
      mean(
        AO_intensity_ratio,
        na.rm = TRUE
      )
    
  )




cluster_profile <-
  
  cluster_profile %>%
  
  mutate(
    
    climate_type = case_when(
      
      
      ENSO + PDO > 0.6 ~
        
        "ENSO-PDO type",
      
      
      AO + NAO > 0.6 ~
        
        "AO-NAO type",
      
      
      TRUE ~
        
        "Multimodal balanced type"
      
    )
    
  )



print(cluster_profile)



# =====================================================
# 7. Add climate classification
# =====================================================


dominance_ratio <-
  
  dominance_ratio %>%
  
  left_join(
    
    cluster_profile %>%
      
      select(
        cluster,
        climate_type
      ),
    
    by="cluster"
    
  )



write.csv(
  
  dominance_ratio,
  
  "climate_dominance_classification.csv",
  
  row.names = FALSE
  
)



write.csv(
  
  cluster_profile,
  
  "climate_type_profiles.csv",
  
  row.names = FALSE
  
)



# =====================================================
# 8. Spatial distribution map
# =====================================================



plot_data <-
  
  dominance_ratio %>%
  
  select(
    
    lon,
    
    lat,
    
    cluster,
    
    climate_type
    
  )




# China boundary

china_province <-
  
  ne_states(
    
    country="China",
    
    returnclass="sf"
    
  )




# Plot

p_map <-
  
  ggplot() +
  
  geom_sf(
    
    data = china_province,
    
    fill = NA,
    
    color = "grey50",
    
    linewidth = 0.4
    
  ) +
  
  
  geom_point(
    
    data = plot_data,
    
    aes(
      
      x = lon,
      
      y = lat,
      
      color = climate_type
      
    ),
    
    size = 4
    
  ) +
  
  
  scale_color_manual(
    
    name = "",
    
    values = c(
      
      "ENSO-PDO type" =
        "#1C3C63",
      
      "AO-NAO type" =
        "#7f280c",
      
      "Multimodal balanced type" =
        "#93C8C0"
      
    )
    
  ) +
  
  
  labs(
    
    x="Longitude",
    
    y="Latitude"
    
  ) +
  
  
  coord_sf(
    
    xlim=c(99,125),
    
    ylim=c(30,46),
    
    expand=FALSE
    
  ) +
  
  
  theme_minimal() +
  
  theme(
    
    panel.grid =
      element_blank(),
    
    
    panel.border =
      element_rect(
        
        color="black",
        
        fill=NA,
        
        linewidth=0.8
        
      ),
    
    
    legend.position =
      "bottom",
    
    
    text =
      element_text(
        size=12
      )
    
  )



print(p_map)



ggsave(
  
  "climate_dominance_map.pdf",
  
  p_map,
  
  width=8,
  
  height=6
  
)



cat(
  
  "Multivariate analysis completed.\n"
  
)