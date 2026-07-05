#Model outputs
#Author: Lizbeth G Amador 


#This script plots brms outputs 

#model type and run # settings 
m = "funct_type" # Could also be nativity and then re-run script
run = "10" # run number out of 10

###########
# Presets
###########
# Load necessary libraries, installing if missing
# note we may not need all these
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}; library(tidyverse)

if (!requireNamespace("brms", quietly = TRUE)) {
  install.packages("brms")
}; library(brms)

if (!requireNamespace("bayesplot", quietly = TRUE)) {
  install.packages("bayesplot")
}; library(bayesplot)

if (!requireNamespace("tidybayes", quietly = TRUE)) {
  install.packages("tidybayes")
}; library(tidybayes)

if (!requireNamespace("marginaleffects", quietly = TRUE)) {
  install.packages("marginaleffects")
}; library(marginaleffects)

if (!requireNamespace("collapse", quietly = TRUE)) {
  install.packages("collapse")
}; library(collapse)

if (!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel")
}; library(parallel)


##############
# Directories
##############

# Note you will need to set this working directory with the GUI: /iplant/home/esokol/esiil_macrophenology/Across_spp_SDM_Scripts/Analysis/Regression/Outputs/Run_10
# graphs <- file.path(getwd(), "Graphs") # writing the graphs will need to be update to go to your local directory

if (!dir.exists(graphs)) {
  dir.create(graphs, recursive = TRUE)
  message("Created directory at: ", graphs)
} else {
  message("Using existing directory at: ", graphs)
}


#########
# Data
#########
# load(file = file.path("/users/PUOM0017/lamad/esiil_macrophenology/Data/L2/Subset_Data_ModelData_NoOutliers.RData"), verbose = TRUE)
load("Subset_Data_ModelData_NoOutliers.RData", verbose = TRUE) # replace this with either training or testing subset data. 
# make sure for the loaded data that it is subsetted further for nativity analysis (hebaceous only)

# Training data
fitb = readRDS(file = paste0("brms_fit_training_Run10_", m, ".rds"))
#Testing data
fitt = readRDS(file = paste0("brms_fit_testing_Run10_", m, ".rds"))

#for file naming 
mod = paste0("Mod", m)

#function to plot both training & testing models
# Fit will be fitt or fitb (training or testing)
# mod will be 'func_type' or 'nativity'
# type will be either 'training' or 'testing'
# run is run number out of ten
fit.plot <- function(fit, mod, type, run){
  message("Plotting!")
  #backing up name
  fit = fit
  #sample for re-predictions 
  n = 10
  
  ####################
  # Model Convergence
  ####################
  pdf(file=file.path(paste0(mod, "_", type, "_R", run, "_brms_Plot_convergence.pdf")))
  #plot dist
  plot(fit)
  # Posterior means + intervals
  mcmc_plot(fit, type = "areas", prob = 0.95)
  dev.off()
  message("Plotted & saved brms convergence plots!")
  
  
  #############
  # Box plots
  #############
  #relationships with pdiff
  #pdiff ~ status 
  #Would be a good idea to make this into a script for nativity and a script for functional type at this point. Note status is not referred to as nativity
  p = plot(conditional_effects(fit, effects = "status"))[[1]]
  p + theme_minimal() +
    theme(plot.background = element_rect(fill = "gray99", color = NA),
      axis.text = element_text(size = 22),
      # axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),
      axis.title = element_text(size = 22),
      plot.title = element_text(size = 22, hjust = 0.5),
      legend.text = element_text(size = 22),
      legend.title = element_text(size = 22))
  
  ggsave(file.path(graphs, paste0(mod, "_", type, "_R", run, "_brms_Plot_pdiff_status_boxplots.png")), last_plot(),
         width = 18, height = 11, units = "in", dpi = 600)
  
  #pdiff ~ functional type
  p = plot(conditional_effects(fit, effects = "functional_type"))[[1]]
  p + theme_minimal() +
    theme(
      plot.background = element_rect(fill = "gray99", color = NA),
      axis.text = element_text(size = 22),
      axis.text.x = element_text(angle = 45, hjust = .9, vjust = .9),
      axis.title = element_text(size = 22),
      plot.title = element_text(size = 22, hjust = 0.5),
      legend.text = element_text(size = 22),
      legend.title = element_text(size = 22)
    )
  
  ggsave(file.path(graphs, paste0(mod, "_", type, "_R", run, "_brms_Plot_pdiff_functype_boxplots.png")), last_plot(),
         width = 16, height = 11, units = "in", dpi = 600)
  
  message("Plotted & saved brms boxplots plots!")
  
  
  ########################### 
  # Plot 1: Coefficient plot 
  ###########################
  #double check font size as it is set really big for a poster
  mcmc_plot(fit, type = "intervals", prob = 0.95, point_size = 5, # make points larger 
            point_color = "black", size = 1.5) + # make interval lines thicker 
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
    theme(text = element_text(size = 22))
  
  ggsave(file.path(graphs, paste0(mod, "_", type, "_R", run, "_brms_Plot_coefficients.png")), last_plot(),
         width = 16, height = 11, units = "in", dpi = 600)
  
  
  #Extract and tidy posterior summaries
  coef_df <- as.data.frame(posterior_summary(fit))
  coef_df$Term <- rownames(coef_df)
  #Filter only fixed effects
  fixef_df <- coef_df %>%
    filter(grepl("^b_", Term)) %>%
    mutate(
      #Term = gsub("^b_", "", Term),
      Sig = ifelse(`Q2.5` > 0 | `Q97.5` < 0, "Significant", "Not Significant"))
  #Plot with ggplot
  ggplot(fixef_df, aes(x = reorder(Term, Estimate), y = Estimate, color = Sig)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
    geom_point(size = 5) +
    geom_errorbar(aes(ymin = `Q2.5`, ymax = `Q97.5`), width = 0.5) +
    coord_flip() +
    labs(x = "", y = "Estimate", title = "") +
    scale_color_manual(values = c("Significant" = "#40B0A6", "Not Significant" = "gray50")) +
    scale_x_discrete(labels = c( # scale_x_discrete labels change axis ticks (variables), may need updating to also add dipseral categories
      "b_Intercept" = "Intercept",
      "b_rdiff" = "Change in Range",
      "b_statusI" = "Status (Invasive)",
      "b_latitude" = "Latitude",
      "b_rdiff:statusI" = "Change in Range × Status (Invasive)",
      "b_rdiff:latitude" = "Change in Range × Latitude",
      "b_functional_typeEvergreenbroadleaf" = "Evergreen Broadleaf",
      "b_functional_typeSemiMevergreenforb" = "Semi-evergreen Forb", 
      "b_functional_typeSemiMevergreenbroadleaf" = "Semi-evergreen Broadleaf",
      "b_functional_typeEvergreenforb" = "Evergreen Forb",
      "b_functional_typeForb" = "Forb",
      "b_functional_typeDroughtdeciduousbroadleaf" = "Drought-deciduous Broadleaf",
      "b_functional_typeSemievergreenforb" = "Semi-evergreen Forb",
      "b_functional_typeSemievergreenbroadleaf" = "Semi-evergreen Broadleaf")) + 
    theme_minimal() +
    theme(plot.background = element_rect(fill = "gray99", color = NA),
          axis.text = element_text(size = 22),
          axis.text.x = element_text(angle = 0, hjust = .4, vjust = 1, size = 22), #est
          axis.text.y = element_text(angle = 0, hjust = .9, vjust = 0.7, size = 22), #terms 
          axis.title = element_text(size = 22),
          plot.title = element_text(size = 1, hjust = 0.5),  # Centered
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 0),
          legend.box.background = element_rect(fill = "white", color = "grey99"),
          legend.position = c(0.75, 0.05),   # inside, nudged left & up
          legend.justification = c("left", "bottom"),
          plot.margin = margin(10, 70, 10, 10))  # add a bit of right margin so legend isn’t clipped
  
  ggsave(file.path(graphs, paste0(mod, "_", type, "_R", run, "_brms_Plot_coefficients.png")), last_plot(),
         width = 16, height = 11, units = "in", dpi = 600, limitsize=FALSE)
  
  rm(coef_df, fixef_df); gc()
  message("Plotted & saved coefficient plot!")
  
  
  ################################# 
  # Plot 2: Simple regression line
  ################################# 
  newdata_grid <- sub.data %>%
    select(rdiff, functional_type, latitude, status, species, dispersal, domain_id) %>%
    distinct() %>%
    slice_sample(n = 1000)
  preds_simple <- add_predicted_draws(fit, newdata = newdata_grid, ndraws = n)
  
  # Plot with posterior predictive intervals
  # add points showing the training or testing data. for nativity you'll need to subset to just have data points from herbaceous species
  model_fit <- preds_simple %>%
    ggplot(aes(x = rdiff, y = .prediction)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50), alpha = 0.5, colour = "black") +
    scale_fill_brewer(palette = "Greys") +
    labs(y = "Change in Phenology\n", x = "\nChange in Range",
      color = "CI",
      fill = "CI") +
    theme_bw() +
    theme(plot.background = element_rect(fill = "gray99", color = NA),
      axis.text = element_text(size = 22),
      axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),
      axis.title = element_text(size = 22),
      plot.title = element_text(size = 22, hjust = 0.5),
      legend.text = element_text(size = 22),
      legend.title = element_text(size = 22),
      legend.position = c(0.92, 0.1),
      legend.box.background = element_rect(fill = "white", color = "grey99"))
  
  ggsave(file.path(graphs, paste0(mod, "_", type, "_R", run, "_brms_Plot_pdiff_rdiff_poster.png")), model_fit,
         width = 16, height = 11, units = "in", dpi = 600)
  
  rm(model_fit, newdata_grid, preds_simple); gc()
  message("Plotted & saved simple regression plot!")
  
  
  ##################################################  
  # Plot 3: Interactions grouped by latitude region 
  ##################################################  
  # Rohit check with Daijiang about if he would like data points added to this one or not.
  sample_data <- sub.data %>%
    mutate(
      lat_region = case_when(
        latitude < 30 ~ "Low (<30°)",
        latitude >= 30 & latitude <= 60 ~ "Mid (30°–60°)",
        latitude > 60 ~ "High (>60°)"))
  sample_data$lat_region <- factor(sample_data$lat_region,
                                   levels = c("Low (<30°)", "Mid (30°–60°)", "High (>60°)"))
  #Add predicted values per lat_region
  preds_lat_region <- sample_data %>%
   group_by(lat_region) %>%
   add_epred_draws(fit, ndraws = n)  #reduce ndraws if needed for speed
  #plot
  ggplot(preds_lat_region, aes(x = rdiff, y = .epred, color = lat_region, fill = lat_region)) +
    stat_lineribbon(.width = c(0.95, 0.8), alpha = 0.25) +
    geom_line(stat = "summary", fun = mean, linewidth = 3) +
    labs(x = "Range Difference (rdiff)", y = "Predicted Phenology Difference (pdiff)",
      color = "Latitude Region",
      fill = "Latitude Region") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "gray99", color = NA),
      axis.text = element_text(size = 22),
      axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),
      axis.title = element_text(size = 22),
      plot.title = element_text(size = 22, hjust = 0.5),
      legend.text = element_text(size = 22),
      legend.title = element_text(size = 22),
      legend.position = c(0.92, 0.1),
      legend.box.background = element_rect(fill = "white", color = "grey99"))
  
  ggsave(file.path(graphs, paste0(mod, "_", type, "_R", run, "_brms_Plot_latitude_interaction.png")), last_plot(),
         width = 16, height = 11, units = "in", dpi = 600)
  
  rm(preds_lat_region); gc()
  message("Plotted & saved latitude interaction plot!")
  
  
  ######################################## 
  # Plot 4: Using marginaleffects package 
  ########################################
  # This plot provides another way to view the interactions
  update_geom_defaults("line", list(linewidth = 3))
  
  #pdiff ~ rdiff*status
  p = plot_predictions(fit, 
                   condition = c("rdiff", "status"),  # x1 on x-axis, lines per x2
  		             ndraws = n)
  p + theme_minimal() +
      theme(plot.background = element_rect(fill = "gray99", color = NA),
        axis.text = element_text(size = 22),
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22, hjust = 0.5),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.position = c(0.92, 0.1),
        legend.box.background = element_rect(fill = "white", color = "grey99"))
  
  ggsave(file.path(graphs, paste0(mod, "_", type, "_R", run, "_brms_Plot_status.png")), last_plot(),
         width = 16, height = 11, units = "in", dpi = 600)
  
  #pdiff ~ rdiff*latitude
  p = plot_predictions(fit, 
                   condition = c("rdiff", "latitude"),  # x1 on x-axis, lines per x2
                   ndraws = n)
  p + theme_minimal() +
    theme(plot.background = element_rect(fill = "gray99", color = NA),
      axis.text = element_text(size = 22),
      axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),
      axis.title = element_text(size = 22),
      plot.title = element_text(size = 22, hjust = 0.5),
      legend.text = element_text(size = 22),
      legend.title = element_text(size = 22),
      legend.position = c(0.92, 0.1),
      legend.box.background = element_rect(fill = "white", color = "grey99"))
  
  ggsave(file.path(graphs, paste0(mod, "_", type, "_R", run, "_brms_Plot_latitude_interaction_cont.png")), last_plot(), 
         width = 16, height = 11, units = "in", dpi = 600)
  
  #pdiff ~ rdiff*functional_type # needs to be updated for nativity
  p = plot_predictions(fit, 
                   condition = c("rdiff", "functional_type"),  # x1 on x-axis, lines per x2
                   ndraws = n)
  p + theme_minimal() +
    theme(plot.background = element_rect(fill = "gray99", color = NA),
      axis.text = element_text(size = 22),
      axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),
      axis.title = element_text(size = 22),
      plot.title = element_text(size = 22, hjust = 0.5),
      legend.text = element_text(size = 22),
      legend.title = element_text(size = 22),
      legend.position = c(0.80, 0.2),
      legend.box.background = element_rect(fill = "white", color = "grey99")) 
  
  ggsave(file.path(graphs, paste0(mod, "_", type, "_R", run, "_brms_Plot_functype.png")), last_plot(),
         width = 16, height = 11, units = "in", dpi = 600)

  message("Plotted & saved marginal error plots!")
  }









