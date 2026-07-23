# =====================================================
# 04_permutation_analysis.R
#
# Quantify the contribution of large-scale climate
# variability modes to regional hydroclimate using
# permutation-based approaches.
#
# Main method:
#   Drop-one conditional permutation
#
# Responses:
#   - precipitation anomaly
#   - temperature anomaly
#   - drought index
#
# =====================================================



# =====================================================
# 1. Load packages
# =====================================================

library(dplyr)
library(tidyr)
library(tibble)
library(VineCopula)



# =====================================================
# 2. Load regional climate surface
# =====================================================


regional_surface <- read.csv(
  "regional_surface.csv",
  stringsAsFactors = FALSE
)



# =====================================================
# 3. Prepare analysis dataset
# =====================================================


climate_perm_df <- regional_surface %>%
  
  select(
    
    PG_s,
    TG_s,
    spei,
    
    ENSO_sim,
    PDO_sim,
    NAO_sim,
    AO_sim,
    SeaIce_sim
    
  ) %>%
  
  drop_na()



# =====================================================
# 4. Define predictor variables
# =====================================================


predictors <- c(
  
  "ENSO_sim",
  "PDO_sim",
  "NAO_sim",
  "AO_sim",
  "SeaIce_sim"
  
)



# =====================================================
# Part I
# Simple permutation method
# (optional comparison)
# =====================================================


permute_lm_effect <- function(
    df,
    response,
    predictors,
    nperm = 1000
){
  
  
  fml <- as.formula(
    
    paste(
      response,
      "~",
      paste(
        predictors,
        collapse = "+"
      )
    )
    
  )
  
  
  lm_base <- lm(
    fml,
    data=df
  )
  
  
  R2_base <-
    
    summary(lm_base)$r.squared
  
  
  
  res <- lapply(
    
    predictors,
    
    function(var){
      
      
      R2_perm <- numeric(nperm)
      
      
      
      for(i in 1:nperm){
        
        
        df_perm <- df
        
        
        df_perm[[var]] <-
          
          sample(
            df_perm[[var]]
          )
        
        
        lm_p <-
          
          lm(
            fml,
            data=df_perm
          )
        
        
        R2_perm[i] <-
          
          summary(lm_p)$r.squared
        
      }
      
      
      
      tibble(
        
        Response=response,
        
        Predictor=var,
        
        R2_orig=R2_base,
        
        R2_perm_mean=
          mean(R2_perm),
        
        Delta_R2=
          R2_base -
          mean(R2_perm),
        
        p_value=
          mean(
            R2_perm >= R2_base
          )
        
      )
      
    }
    
  )
  
  
  bind_rows(res)
  
}




# =====================================================
# Part II
# Main analysis:
# Drop-one conditional permutation
# =====================================================



drop_one_conditional_perm <- function(
    df,
    response,
    predictors,
    nperm = 1000
){
  
  
  # Full model
  
  fml_full <-
    
    as.formula(
      
      paste(
        response,
        "~",
        paste(
          predictors,
          collapse="+"
        )
      )
      
    )
  
  
  lm_full <-
    
    lm(
      fml_full,
      data=df
    )
  
  
  R2_full <-
    
    summary(lm_full)$r.squared
  
  
  
  
  res <- lapply(
    
    predictors,
    
    function(var){
      
      
      
      # remove one predictor
      
      preds_drop <-
        
        setdiff(
          predictors,
          var
        )
      
      
      
      fml_drop <-
        
        as.formula(
          
          paste(
            response,
            "~",
            paste(
              preds_drop,
              collapse="+"
            )
          )
          
        )
      
      
      
      lm_drop <-
        
        lm(
          fml_drop,
          data=df
        )
      
      
      
      R2_drop <-
        
        summary(lm_drop)$r.squared
      
      
      
      Delta_R2 <-
        
        R2_full -
        R2_drop
      
      
      
      
      # Conditional permutation
      # Preserve correlations among predictors
      
      
      lm_var <-
        
        lm(
          
          df[[var]] ~
            
            .,
          
          data =
            df[
              ,
              preds_drop,
              drop=FALSE
            ]
          
        )
      
      
      
      fitted_var <-
        
        fitted(lm_var)
      
      
      residual_var <-
        
        residuals(lm_var)
      
      
      
      Delta_perm <-
        
        numeric(nperm)
      
      
      
      for(i in 1:nperm){
        
        
        
        df_perm <- df
        
        
        
        df_perm[[var]] <-
          
          fitted_var +
          sample(
            residual_var
          )
        
        
        
        lm_full_perm <-
          
          lm(
            fml_full,
            data=df_perm
          )
        
        
        
        lm_drop_perm <-
          
          lm(
            fml_drop,
            data=df_perm
          )
        
        
        
        Delta_perm[i] <-
          
          summary(lm_full_perm)$r.squared -
          summary(lm_drop_perm)$r.squared
        
        
      }
      
      
      
      p_value <-
        
        (
          sum(
            Delta_perm >= Delta_R2
          ) + 1
        ) /
        (nperm + 1)
      
      
      
      tibble(
        
        Response=response,
        
        Predictor=var,
        
        R2_full=R2_full,
        
        R2_drop=R2_drop,
        
        Delta_R2=Delta_R2,
        
        p_value=p_value
        
      )
      
    }
    
  )
  
  
  
  bind_rows(res)
  
}





# =====================================================
# 5. Run conditional permutation analysis
# =====================================================



perm_PG <-
  
  drop_one_conditional_perm(
    
    climate_perm_df,
    
    response="PG_s",
    
    predictors=predictors,
    
    nperm=1000
    
  )



perm_TG <-
  
  drop_one_conditional_perm(
    
    climate_perm_df,
    
    response="TG_s",
    
    predictors=predictors,
    
    nperm=1000
    
  )



perm_SPEI <-
  
  drop_one_conditional_perm(
    
    climate_perm_df,
    
    response="spei",
    
    predictors=predictors,
    
    nperm=1000
    
  )




perm_climate_results <-
  
  bind_rows(
    
    perm_PG,
    
    perm_TG,
    
    perm_SPEI
    
  )



print(
  perm_climate_results,
  n=Inf,
  width=Inf
)



write.csv(
  
  perm_climate_results,
  
  "LCVM_regional_climate_contribution.csv",
  
  row.names=FALSE
  
)



# =====================================================
# Part III
# Dependency assessment among LCVM predictors
# =====================================================



predictor_data <-
  
  climate_perm_df[, predictors]



# Pearson

pearson_matrix <-
  
  cor(
    
    predictor_data,
    
    method="pearson",
    
    use="pairwise.complete.obs"
    
  )



# Spearman

spearman_matrix <-
  
  cor(
    
    predictor_data,
    
    method="spearman",
    
    use="pairwise.complete.obs"
    
  )



# Kendall

kendall_matrix <-
  
  cor(
    
    predictor_data,
    
    method="kendall",
    
    use="pairwise.complete.obs"
    
  )



write.csv(
  pearson_matrix,
  "LCVM_Pearson_correlation.csv"
)


write.csv(
  spearman_matrix,
  "LCVM_Spearman_correlation.csv"
)


write.csv(
  kendall_matrix,
  "LCVM_Kendall_correlation.csv"
)



# =====================================================
# Part IV
# Vine Copula dependence bootstrap
# =====================================================


U_large <-
  
  climate_perm_df %>%
  
  select(
    all_of(predictors)
  ) %>%
  
  mutate(
    
    across(
      everything(),
      ~rank(.)/(length(.)+1)
    )
    
  )



vine_large <-
  
  RVineStructureSelect(
    
    as.matrix(U_large),
    
    type="CVine"
    
  )



tau_matrix <-
  
  RVinePar2Tau(
    vine_large
  )



write.csv(
  
  tau_matrix,
  
  "LCVM_VineCopula_Kendall_tau.csv"
  
)



cat(
  "Permutation analysis completed.\n"
)
# =====================================================
# 5. Convert permutation results into dominance ratio
# =====================================================


# 选择用于分类的响应变量
# 例如区域气候综合响应 SPEI

dominance_input <- perm_results %>%
  
  filter(
    Response == "spei"
  )


# 将 Delta_R2 转换为相对贡献比例

dominance_ratio <- dominance_input %>%
  
  group_by(Response) %>%
  
  mutate(
    
    ratio =
      Delta_R2 /
      sum(Delta_R2, na.rm = TRUE)
    
  ) %>%
  
  ungroup()



# 转换为宽格式

dominance_ratio <- dominance_ratio %>%
  
  select(
    Predictor,
    ratio
  ) %>%
  
  tidyr::pivot_wider(
    
    names_from = Predictor,
    
    values_from = ratio
    
  )



# 修改变量名，匹配05脚本

dominance_ratio <- dominance_ratio %>%
  
  rename(
    
    ENSO_intensity_ratio =
      ENSO_sim,
    
    PDO_intensity_ratio =
      PDO_sim,
    
    NAO_intensity_ratio =
      NAO_sim,
    
    AO_intensity_ratio =
      AO_sim,
    
    SeaIce_intensity_ratio =
      SeaIce_sim
    
  )



write.csv(
  
  dominance_ratio,
  
  "dominance_ratio.csv",
  
  row.names = FALSE
  
)