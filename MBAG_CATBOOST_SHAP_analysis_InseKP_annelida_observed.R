### ESB CATBoost machine learning model 
### Marginal effect of stratifiers on prediction of target variables
### Determines Relative importance of predictors
### Creates beeswarm plots & dependence plots for analysis and exploration
### Makes predictive model based on training dataset
### Evaluation of predictive capacity of model using PQ metrics (R², bias, RMSE)
### Bruno De Vos 2025-09-17

library(dplyr)
library(ggplot2)
library(catboost)
library(shapviz)


## SHAP = SHAP stands for SHapley Additive exPlanations.
# --- shapviz integration ---
# Wrapper for using shapviz with CatBoost
shapviz.catboost.Model <- function(object, X_pred, X = X_pred, collapse = NULL, ...) {
  if (!inherits(X_pred, "catboost.Pool")) {
    X_pred <- catboost.load_pool(X_pred)
  }
  S <- catboost.get_feature_importance(object, X_pred, type = "ShapValues", ...)
  pp <- ncol(S)
  baseline <- S[1, pp]
  S <- S[, -pp, drop = FALSE]
  colnames(S) <- colnames(X)
  shapviz(S, X = X, baseline = baseline, collapse = collapse)
}



source(file = "C:/R_scripts/R_PQ/PQD_function.R")

# Set working directory ---PUT YOUR OWN WD here --- 
setwd("C:/data/MBAG/")


#Load soil data  ------
DS_MBAG <- read.csv("./mbag_combined_dataframe_metadata_wide_nematoda_annelida.csv")
names(DS_MBAG)
str(DS_MBAG)

## target variable (absolute)
DS_MBAG$InseKP_annelida_observed
## or target value log-transformed (to stabilize variance and make the distribution more Gaussian-like)
# when you want the model to focus on relative rather than absolute differences.


## dataset with non NA DS_MBAG_Anne
DS_MBAG_Anne <- DS_MBAG[!is.na(DS_MBAG$InseKP_annelida_observed), ]
summary(DS_MBAG_Anne)

# attributes
names(DS_MBAG_Anne)

### ###
# training dataset 
### ###

# Prepare the dataset
## target response variabele 
Anne_Richness<-as.numeric(DS_MBAG_Anne$InseKP_annelida_observed) 
summary(Anne_Richness)
hist(Anne_Richness)
Anne_Richness_lt<-log(Anne_Richness+1)
hist(Anne_Richness_lt)


## recode
unique(DS_MBAG_Anne$Diepte)
DS_MBAG_Anne$Diepte<-ifelse(DS_MBAG_Anne$Diepte==45960,"10-30", DS_MBAG_Anne$Diepte)


# predictors = features 
# VERY important: categorical predictors as factors, continuous variables as numeric 
s1_feat <- data.frame(
  LogTotCount =as.numeric(log(DS_MBAG_Anne$InseKP_annelida_total_count+1)),      # with total count "correction"
  CD = as.numeric(DS_MBAG_Anne$Cdensity),
  ND = as.numeric(DS_MBAG_Anne$Ndensity),
  pH = as.numeric(DS_MBAG_Anne$pH_KCl),
  BD = as.numeric(DS_MBAG_Anne$BD),
  TOC = as.numeric(DS_MBAG_Anne$TOC),  
  TN = as.numeric(DS_MBAG_Anne$TN), 
  CNR = as.numeric(DS_MBAG_Anne$C_N_stockbased),
  SWC = as.numeric(DS_MBAG_Anne$SWCvol),
  CLAY = as.numeric(DS_MBAG_Anne$Textuur_kleifractie),
  LU = as.factor(DS_MBAG_Anne$Landgebruik_MBAG),
  DEPTH = as.factor(DS_MBAG_Anne$Diepte),
  TEX = as.factor(DS_MBAG_Anne$Textuurklasse),
  DRAIN = as.factor(DS_MBAG_Anne$Drainageklasse)
)

## summary predictors
summary(s1_feat)

## check continuous / categorical predictors
str(s1_feat)

#View(s1_feat)

# Create a data matrix for CatBoost
s1_pool <- catboost.load_pool(data = s1_feat, label = Anne_Richness)


# Important to avoid overfitting !!
# Set parameters with early stopping using cross-validation 

params <- list(
  loss_function = "RMSE",
  iterations = 5000,              # Start with a high number
  od_type = "Iter",               # Use iteration-based early stopping
  od_wait = 30,                   # Wait 30 iterations before stopping
  allow_writing_files = FALSE,
  random_seed = 123                # For reproducibility
)

# Perform cross-validation to find optimal number of iterations
cv_results <- catboost.cv(
  pool = s1_pool,
  params = params,
  fold_count = 5,                 # 5-fold cross-validation
  type = "Classical",             # Classical Cross validation 
  partition_random_seed = 123
)

# Extract best iteration based on minimum RMSE
best_iter <- which.min(cv_results$test.RMSE.mean)
cat("Best number of iterations based on CV:", best_iter, "\n")

# Train final model using best teration
params$iterations <- best_iter
Anne_Richness_CBmodel  <- catboost.train(s1_pool, params = params)
summary(Anne_Richness_CBmodel)    ## this is your CATBOOST model for training dataset


# Visualize with shapviz
# Make SHAPVIZ object
Anne_shap <- shapviz(Anne_Richness_CBmodel, X_pred = s1_pool, X = s1_feat)
saveRDS(Anne_shap, file = "./Anne_shap_model_Ann.rds")   ## save your SHAP object for later SHAP analysis

#read it with Anne_shap <- readRDS("Anne_shap_model.rds")
# mean SHAP values for each predictor (relative importance)
Ann_imp<-sv_importance(Anne_shap, kind="no")    # Table with average Anne_shap values for each feature
#write.csv2(Ann_imp,"./CB_Ann_importance.csv", row.names = TRUE)
round(Ann_imp,2)

#bar graph
sv_importance(Anne_shap, kind = "bar")

# check potential interactions between predictors 
# This measures how much variability in the SHAP values of v is explained by ⁠v'⁠, after accounting for v.
LogTotCount_inter<-potential_interactions(Anne_shap,"LogTotCount",adjusted=TRUE) ; LogTotCount_inter
DEPTH_inter<-potential_interactions(Anne_shap,"DEPTH") ; DEPTH_inter
LU_inter<-potential_interactions(Anne_shap,"LU") ; LU_inter

# overall importance map
sv_importance(Anne_shap, kind = "beeswarm", show_numbers = TRUE)


## dependence plots of interest

#TotCount
summary(s1_feat$LogTotCount)

### dependance 
p <- sv_dependence(Anne_shap, "LogTotCount", color_var = "pH")
p +
  scale_color_gradientn(
    colors = c("darkblue", "green","red"),
    limits = c(2, 10)  # Set scale 
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1.2) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = TRUE,                 # show confidence interval
    color = "blue",            # line color
    fill = "grey70",           # ribbon fill color
    alpha = 0.4,               # transparency of ribbon
    linewidth = 1.2
  ) +
  ggtitle("Annelida richness SHAP-based Dependence Plot on log total_count with pH scaling") +
  scale_x_continuous(limits = c(0, 16)) +  # <-- set your desired x-axis range
  labs(
    x = "Annelids_log(total_count+1)",   # <-- custom X axis label
    y = "SHAP value"                        # <-- custom Y axis label
  ) +
  theme_minimal()




#sv_dependence(Anne_shap, colnames(LI_feat)) ## all features
# C:N ratio
p <- sv_dependence(Anne_shap, "CNR", color_var = "pH")
p +
  scale_color_gradientn(
    colors = c("darkblue", "green","red"),
    limits = c(2, 10)  # Set scale 
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1.2) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = TRUE,                 # show confidence interval
    color = "blue",            # line color
    fill = "grey70",           # ribbon fill color
    alpha = 0.4,               # transparency of ribbon
    linewidth = 1.2
  ) +
  ggtitle("Annelida richness SHAP-based Dependence Plot on C:N ratio with pH scaling") +
  scale_x_continuous(limits = c(5, 30)) +  # <-- set your desired x-axis range
  labs(
    x = "C:N ration)",   # <-- custom X axis label
    y = "SHAP value"                        # <-- custom Y axis label
  ) +
  theme_minimal()


# pH
p <- sv_dependence(Anne_shap, "pH", color_var = "CNR")
p +
  scale_color_gradientn(
    colors = c("darkblue", "green","red"),
    limits = c(0, 30)  # Set scale 
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1.2) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = TRUE,                 # show confidence interval
    color = "blue",            # line color
    fill = "grey70",           # ribbon fill color
    alpha = 0.4,               # transparency of ribbon
    linewidth = 1.2
  ) +
  ggtitle("Annelida richness SHAP-based Dependence Plot on pH with C:N ratio scaling") +
  scale_x_continuous(limits = c(2.5, 8)) +  # <-- set your desired x-axis range
  labs(
    x = "pH value",   # <-- custom X axis label
    y = "SHAP value"                        # <-- custom Y axis label
  ) +
  theme_minimal()


# CATEGORICAL DEPTH
p <- sv_dependence(Anne_shap, "DEPTH", color_var = "LU")
p + geom_hline(
  yintercept = 0,
  color = "black",
  linetype = "dashed",
  linewidth = 1.2
) +
  ggtitle("InseKP_annelida_observed richness SHAP-based Dependence Plot for landuse with Depth colours") +
  labs(
    x = "Landuse category",   # <-- custom X axis label
    y = "SHAP value"                        # <-- custom Y axis label
  ) + theme_minimal()




# CATEGORICAL LU
p <- sv_dependence(Anne_shap, "LU", color_var = "DEPTH")
p + geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  ggtitle("InseKP_annelida_observed richness SHAP-based Dependence Plot for landuse with Depth colours") +
  labs(
    x = "Landuse category",   # <-- custom X axis label
    y = "SHAP value"                        # <-- custom Y axis label
  ) + theme_minimal()




# Carbon density
p <- sv_dependence(Anne_shap, "CD", color_var = "pH")
p +
  scale_color_gradientn(
    colors = c("darkblue", "green","red"),
    limits = c(2, 10)  # Set scale 
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1.2) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = TRUE,                 # show confidence interval
    color = "blue",            # line color
    fill = "grey70",           # ribbon fill color
    alpha = 0.4,               # transparency of ribbon
    linewidth = 1.2
  ) +
  ggtitle("Annelida richness SHAP-based Dependence Plot on Carbon Density with pH scaling") +
  scale_x_continuous(limits = c(0, 8)) +  # <-- set your desired x-axis range
  labs(
    x = "Carbon Density (t·C ha⁻¹ cm⁻¹)",   # <-- custom X axis label
    y = "SHAP value"                        # <-- custom Y axis label
  ) +
  theme_minimal()




# Total Nitrogen
summary(s1_feat$TN)
# plotting
p <- sv_dependence(Anne_shap, "TN", color_var = "CD")
p +
  scale_color_gradientn(
    colors = c("darkblue", "green","red"),
    limits = c(0, 10)  # Set scale 
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1.2) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = TRUE,                 # show confidence interval
    color = "blue",            # line color
    fill = "grey70",           # ribbon fill color
    alpha = 0.4,               # transparency of ribbon
    linewidth = 1.2
  ) +
  ggtitle("Annelida richness SHAP-based Dependence Plot on Total Nitrogen with CD scaling") +
  scale_x_continuous(limits = c(0, 0.8)) +  # <-- set your desired x-axis range
  labs(
    x = "Total Nitrogen (%)",   # <-- custom X axis label
    y = "SHAP value"                        # <-- custom Y axis label
  ) +
  theme_minimal()



# pH
summary(s1_feat$pH)
# plotting
p <- sv_dependence(Anne_shap, "pH", color_var = "LU" )
p +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1.2) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = TRUE,                 # show confidence interval
    color = "blue",            # line color
    fill = "grey70",           # ribbon fill color
    alpha = 0.4,               # transparency of ribbon
    linewidth = 1.2
  ) +
  ggtitle("Annelida richness SHAP-based Dependence Plot on Soil pH with CD scaling") +
  scale_x_continuous(limits = c(2.8, 8.2)) +  # <-- set your desired x-axis range
  labs(
    x = "Soil pH-KCl",   # <-- custom X axis label
    y = "SHAP value"                        # <-- custom Y axis label
  ) +
  theme_minimal()


#  TOC
summary(s1_feat$TOC)
# plotting
p <- sv_dependence(Anne_shap, "TOC", color_var = "LU" )
p +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1.2) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = TRUE,                 # show confidence interval
    color = "blue",            # line color
    fill = "grey70",           # ribbon fill color
    alpha = 0.4,               # transparency of ribbon
    linewidth = 1.2
  ) +
  ggtitle("Annelida richness SHAP-based Dependence Plot on Soil TOC with LU scaling") +
  scale_x_continuous(limits = c(0.1, 12)) +  # <-- set your desired x-axis range
  labs(
    x = "TOC   (%)",   # <-- custom X axis label
    y = "SHAP value"                        # <-- custom Y axis label
  ) +
  theme_minimal()



#####################
## waterfall plots ##
#####################

# Check predict specific plots
sv_waterfall(Anne_shap, row_id=7) 
# for landuses
par(mfrow=c(3,2))

unique(Anne_shap$X$LU)

p<-sv_waterfall(Anne_shap, row_id = Anne_shap$X$LU=="Akker")
p+ggtitle("Waterfall plot Akker (n=178)")+
theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

p<-sv_waterfall(Anne_shap, row_id = Anne_shap$X$LU=="Tijdelijk grasland")
p+ggtitle("Waterfall plot Tijdelijk grasland (n=70)")+
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

p<-sv_waterfall(Anne_shap, row_id = Anne_shap$X$LU=="Blijvend grasland")
p+ggtitle("Waterfall plot Blijvend grasland (n=95)")+
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

p<-sv_waterfall(Anne_shap, row_id = Anne_shap$X$LU=="Residentieel grasland")
p+ggtitle("Waterfall plot Residentieel grasland (n=127)")+
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

p<-sv_waterfall(Anne_shap, row_id = Anne_shap$X$LU=="Natuurgrasland")
p+ggtitle("Waterfall plot Natuurgrasland  (n=92)")+
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))


# for topsoil
p<-sv_waterfall(Anne_shap, row_id = Anne_shap$X$DEPTH=="0-10")
p+ggtitle("Waterfall plot for 0-10 cm layer")+
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

# for 10-30 cm layer
p<-sv_waterfall(Anne_shap, row_id = Anne_shap$X$DEPTH=="10-30")
p+ggtitle("Waterfall plot for 10-30 cm layer")+
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))


#row_id links to specific row (plot_id)
# force plot allows to average over rows (plot_ids)
Anne_shap
sv_force(Anne_shap, row_id=1:586)   ### aggregate over all plots 1:


######################################################################################################
### predictive capacity of CATBOOST trained model

p_S1_Ann_train <- catboost.predict(Anne_Richness_CBmodel, s1_pool);length(p_S1_Ann_train)
length(p_S1_Ann_train)

par(mfrow=c(1,1))

# Explained variance by CATBOOST model
# Calculate R-squared (variance explained) of trained model
Observed_Richness<-Anne_Richness    ## observed value
Predicted_Richness<-p_S1_Ann_train  ## predicted value

length(Observed_Richness);length(Predicted_Richness)

SS_res <- sum((Observed_Richness - Predicted_Richness)^2)
SS_tot <- sum((Observed_Richness - mean(Observed_Richness))^2)
R_squared <- 1 - (SS_res / SS_tot) 
cat("AWC CATBOOST Total Variance Explained (R-squared):", round(R_squared, 2), "\n")

# prediction quality
PQD(Predicted_Richness,Observed_Richness,2)

range(Observed_Richness)
range(Predicted_Richness)


tmax<-max(c(Observed_Richness,Predicted_Richness))
plot(Observed_Richness,Predicted_Richness, pch=16, xlim=c(0,tmax), ylim=c(0,tmax),
     main="CATBOOST predictions of Annelid Richness on MBAG plots")
abline(0,1, lty=2, lwd=2,col="blue")
lines(lowess(Observed_Richness,Predicted_Richness), col="red", lwd=3)
legend("bottomright", legend=c("1:1 line","Lowess curve"),lwd=c(2,3), col=c("blue","red"))


## summary stats 
length(Observed_Richness)   ### number of observed values
summary(Observed_Richness)  ### observed/measured values 
summary(Predicted_Richness)  ### predicted values by model

### from here you can analyse the residuals (p-o) for
# correlation/covariate analysis 
# geostatistical (spatial) analysis (kriging, ...)



