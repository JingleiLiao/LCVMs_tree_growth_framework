# =====================================================
# 03_LSTM_climate_growth_modeling.R
#
# Build long short-term memory (LSTM) models
# linking climate variability modes and tree growth.
#
# Input:
#   regional_surface.csv
#
# Output:
#   lstm_model.keras
#   lstm_performance.csv
#   lstm_predictions.csv
#
# =====================================================



# =====================================================
# 1. Load packages
# =====================================================

library(dplyr)
library(tidyr)
library(abind)
library(keras)
library(tensorflow)



# =====================================================
# 2. Load prepared regional climate-growth dataset
# =====================================================


regional_surface <- read.csv(
  "regional_surface.csv",
  stringsAsFactors = FALSE
)



# =====================================================
# 3. Prepare LSTM dataset
# =====================================================


tree_surface_lstm <- regional_surface %>%
  
  arrange(
    lon,
    lat,
    year
  ) %>%
  
  group_by(
    lon,
    lat
  ) %>%
  
  mutate(
    
    PG_z =
      as.numeric(scale(PG)),
    
    
    TG_z =
      as.numeric(scale(TG)),
    
    
    SPEIG_z =
      as.numeric(scale(SPEIG)),
    
    
    ENSO_z =
      as.numeric(scale(ENSO_sim)),
    
    
    PDO_z =
      as.numeric(scale(PDO_sim)),
    
    
    NAO_z =
      as.numeric(scale(NAO_sim)),
    
    
    AO_z =
      as.numeric(scale(AO_sim)),
    
    
    SeaIce_z =
      as.numeric(scale(SeaIce_sim)),
    
    
    RWI_z =
      as.numeric(scale(rwi_obs))
    
  ) %>%
  
  ungroup()



# =====================================================
# 4. Define LSTM variables
# =====================================================


input_vars <- c(
  
  "PG_z",
  "TG_z",
  "SPEIG_z",
  
  "ENSO_z",
  "PDO_z",
  "NAO_z",
  "AO_z",
  "SeaIce_z"
  
)


output_var <- "RWI_z"


lag <- 6



# =====================================================
# 5. Generate LSTM sequences
# =====================================================


data_list <- tree_surface_lstm %>%
  
  arrange(
    lon,
    lat,
    year
  ) %>%
  
  group_split(
    lon,
    lat
  )



make_lstm_array <- function(
    df,
    input_vars,
    output_var,
    lag
){
  
  df <- df %>%
    drop_na(
      all_of(
        c(
          input_vars,
          output_var
        )
      )
    )
  
  
  if(
    nrow(df) <= lag
  ){
    return(NULL)
  }
  
  
  X <- array(
    NA,
    dim=c(
      nrow(df)-lag,
      lag,
      length(input_vars)
    )
  )
  
  
  Y <- numeric(
    nrow(df)-lag
  )
  
  
  for(i in 1:(nrow(df)-lag)){
    
    X[i,,] <-
      
      as.matrix(
        df[
          i:(i+lag-1),
          input_vars
        ]
      )
    
    
    Y[i] <-
      
      df[[output_var]][i+lag]
    
  }
  
  
  list(
    X=X,
    Y=Y
  )
  
}



X_all <- NULL

Y_all <- NULL



for(df in data_list){
  
  tmp <-
    
    make_lstm_array(
      df,
      input_vars,
      output_var,
      lag
    )
  
  
  if(is.null(tmp)){
    next
  }
  
  
  X_all <-
    
    abind(
      X_all,
      tmp$X,
      along=1
    )
  
  
  Y_all <-
    
    c(
      Y_all,
      tmp$Y
    )
  
}



# =====================================================
# 6. Construct LSTM model
# =====================================================


set.seed(123)

tensorflow::tf$random$set_seed(123)



model <- keras_model_sequential() %>%
  
  layer_lstm(
    units = 32,
    input_shape =
      c(
        lag,
        length(input_vars)
      )
  ) %>%
  
  layer_dense(
    units = 1
  )



model %>%
  
  compile(
    
    optimizer =
      optimizer_adam(
        learning_rate = 0.001
      ),
    
    loss = "mse"
    
  )



# =====================================================
# 7. Train model
# =====================================================


history <- model %>%
  
  fit(
    
    X_all,
    Y_all,
    
    epochs = 80,
    
    batch_size = 64,
    
    validation_split = 0.2,
    
    verbose = 1
    
  )



# =====================================================
# 8. Model evaluation
# =====================================================


prediction <-
  
  as.numeric(
    model %>% predict(X_all)
  )


valid_idx <-
  
  complete.cases(
    Y_all,
    prediction
  )


Y_valid <-
  
  Y_all[valid_idx]


prediction_valid <-
  
  prediction[valid_idx]



lstm_performance <-
  
  data.frame(
    
    R2 =
      cor(
        Y_valid,
        prediction_valid
      )^2,
    
    
    RMSE =
      sqrt(
        mean(
          (prediction_valid-Y_valid)^2
        )
      ),
    
    
    Bias =
      mean(
        prediction_valid-Y_valid
      )
    
  )



print(
  lstm_performance
)



# =====================================================
# 9. Save outputs
# =====================================================


save_model(
  
  model,
  
  "lstm_model.keras"
  
)



write.csv(
  
  lstm_performance,
  
  "lstm_performance.csv",
  
  row.names = FALSE
  
)



write.csv(
  
  data.frame(
    
    observed = Y_valid,
    
    predicted = prediction_valid
    
  ),
  
  "lstm_predictions.csv",
  
  row.names = FALSE
  
)



cat(
  "LSTM modeling completed.\n"
)