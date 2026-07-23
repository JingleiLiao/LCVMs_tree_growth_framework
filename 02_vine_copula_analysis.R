# =====================================================
# 02_vine_copula_analysis.R
#
# Construct nonlinear dependence structures among
# large-scale climate variability modes (LCVMs)
# using Vine Copula models and generate coupled
# climate variability surfaces.
#
# Input:
#   - STD : tree-ring chronology with coordinates
#   - AI  : large-scale climate variability indices
#   - T   : regional climate variables
#
# Output:
#   - vine_large
#   - sim_large
#   - regional_surface
#
# =====================================================


# =====================================================
# 1. Load packages
# =====================================================

library(dplyr)
library(tidyr)
library(VineCopula)



# =====================================================
# 2. Load prepared datasets
# =====================================================

# Users should modify paths according to local systems

STD <- read.csv(
  "tree_ring_chronology.csv",
  stringsAsFactors = FALSE
)


AI <- read.csv(
  "large_scale_climate_modes.csv",
  stringsAsFactors = FALSE
)


T <- read.csv(
  "regional_climate.csv",
  stringsAsFactors = FALSE
)



# =====================================================
# 3. Define climate mode intensity function
# =====================================================

mode_intensity <- function(df, vars){
  
  scale(
    rowMeans(
      dplyr::select(df, all_of(vars)),
      na.rm = TRUE
    )
  )[,1]
  
}



# =====================================================
# 4. Construct spatial-temporal climate variability
#    surface
# =====================================================


years <- STD %>%
  distinct(year)


locations <- STD %>%
  distinct(lon, lat)


grid <- tidyr::crossing(
  years,
  locations
)



large_scale_surface <- grid %>%
  
  left_join(
    AI,
    by = "year"
  ) %>%
  
  mutate(
    lon = lon.x,
    lat = lat.x
  ) %>%
  
  select(
    -lon.x,
    -lat.x,
    -lon.y,
    -lat.y
  )



# =====================================================
# 5. Calculate LCVM intensity indices
# =====================================================


large_scale_surface <- large_scale_surface %>%
  
  mutate(
    
    ENSO_intensity =
      scale(
        rowMeans(
          select(.,
                 ENSOMIEV29:ENSOMIEV32,
                 ensomiev21:ensomiev29),
          na.rm = TRUE
        )
      )[,1],
    
    
    PDO_intensity =
      mode_intensity(
        .,
        c(
          "PDO9","PDO10",
          "PDO11","PDO12",
          "pdo1","pdo2",
          "pdo3","pdo4",
          "pdo5","pdo6",
          "pdo7","pdo8",
          "pdo9","pdoG"
        )
      ),
    
    
    NAO_intensity =
      mode_intensity(
        .,
        c(
          "NAO9","NAO10",
          "NAO11","NAO12",
          "Nao1","Nao2",
          "Nao3","Nao4",
          "Nao5","Nao6",
          "Nao7","Nao8",
          "Nao9","GNAO"
        )
      ),
    
    
    AO_intensity =
      mode_intensity(
        .,
        c(
          "AO9","AO10",
          "AO11","AO12",
          "Ao1","Ao2",
          "Ao3","Ao4",
          "Ao5","Ao6",
          "Ao7","Ao8",
          "Ao9","GAO"
        )
      ),
    
    
    SeaIce_intensity =
      mode_intensity(
        .,
        c(
          "HAIBING9",
          "HAIBING10",
          "HAIBING11",
          "HAIBING12",
          "haibing1",
          "haibing2",
          "haibing3",
          "haibing4",
          "haibing5",
          "haibing6",
          "haibing7",
          "haibing8",
          "haibing9",
          "Ghaibing"
        )
      )
    
  )



# =====================================================
# 6. Combine regional climate and tree growth data
# =====================================================


full_surface <- large_scale_surface %>%
  
  left_join(
    T,
    by = c(
      "year",
      "lon",
      "lat"
    )
  ) %>%
  
  left_join(
    STD,
    by = c(
      "year",
      "lon",
      "lat"
    )
  ) %>%
  
  filter(
    !is.na(PG),
    !is.na(TG),
    !is.na(SPEIG)
  ) %>%
  
  mutate(
    
    PG_s =
      as.numeric(scale(PG)),
    
    
    TG_s =
      as.numeric(scale(TG)),
    
    
    SPEIG_s =
      as.numeric(scale(SPEIG))
    
  )



# =====================================================
# 7. Transform LCVM variables into copula space
# =====================================================


U_large <- large_scale_surface %>%
  
  select(
    ENSO_intensity,
    PDO_intensity,
    NAO_intensity,
    AO_intensity,
    SeaIce_intensity
  ) %>%
  
  mutate(
    across(
      everything(),
      ~ rank(.) /
        (length(.) + 1)
    )
  ) %>%
  
  na.omit()



# =====================================================
# 8. Construct Vine Copula model
# =====================================================


vine_large <- RVineStructureSelect(
  as.matrix(U_large),
  type = "CVine"
)



# =====================================================
# 9. Generate coupled LCVM simulations
# =====================================================


set.seed(123)


n_sim <- nrow(full_surface)


sim_large <- RVineSim(
  n_sim,
  vine_large
) %>%
  as.data.frame()


colnames(sim_large) <- colnames(U_large)



# Add spatial-temporal information

sim_large <- cbind(
  
  sim_large,
  
  full_surface %>%
    select(
      year,
      lon,
      lat
    )
  
)



# Rename simulated LCVM variables

sim_large <- sim_large %>%
  
  rename(
    
    ENSO_sim =
      ENSO_intensity,
    
    PDO_sim =
      PDO_intensity,
    
    NAO_sim =
      NAO_intensity,
    
    AO_sim =
      AO_intensity,
    
    SeaIce_sim =
      SeaIce_intensity
    
  ) %>%
  
  distinct(
    year,
    lon,
    lat,
    .keep_all = TRUE
  )



# =====================================================
# 10. Construct regional climate surfaces
# =====================================================


w <- rep(
  1/5,
  5
)


names(w) <- c(
  "ENSO_sim",
  "PDO_sim",
  "NAO_sim",
  "AO_sim",
  "SeaIce_sim"
)



regional_surface <- full_surface %>%
  
  left_join(
    sim_large,
    by = c(
      "year",
      "lon",
      "lat"
    )
  ) %>%
  
  mutate(
    
    precip_anom =
      PG_s +
      as.matrix(
        select(., names(w))
      ) %*% w +
      rnorm(
        n(),
        0,
        0.3
      ),
    
    
    temp_anom =
      TG_s +
      as.matrix(
        select(., names(w))
      ) %*% w +
      rnorm(
        n(),
        0,
        0.3
      ),
    
    
    spei =
      scale(
        precip_anom -
          0.4 * temp_anom +
          SPEIG_s
      )[,1]
    
  )



# =====================================================
# 11. Save outputs
# =====================================================


write.csv(
  regional_surface,
  "regional_surface.csv",
  row.names = FALSE
)


save(
  vine_large,
  file = "vine_large.RData"
)



# =====================================================
# 12. Check output
# =====================================================


head(
  regional_surface %>%
    select(
      year,
      lon,
      lat,
      precip_anom,
      temp_anom,
      spei
    )
)


cat(
  "Vine Copula analysis completed.\n"
)