## BRMS linear regression with uncertainty
## Authors: Tong Qiu, Lizbeth G Amador


#Libraries
if(!require("xfun")) install.packages("xfun")
xfun::pkg_attach2(c("tidyverse", "brms", "caret", "dggridR"))

#########
# The models were run 10 times to ensure that the spatially stratified random selection of points
# did not generate any bias in the model inferences. This spatial subsetting was done 
# to alleviate computational burdens when running these models.
run = 1 # This is manually updated as runs 1-10 of the model are implemented.
#########

#######
# Data
#######


#Model data
load(file = "data-store/home/shared/esiil/macrophenology/Across_spp_SDM/Data/L2/Subset_ModelData_NoOutliers.RData", verbose = TRUE) # variables of importance: pdiff, psd, rdiff, rsd, domain_id, funct_tpe, nativity, dispersal
#NEON Domain adjaceny table -- Important to read as a matrix!
# cat(1:4, file = "test.txt")
dom.mat <- as.matrix(read.csv("data-store/home/shared/esiil/macrophenology/Across_spp_SDM/Data/L2/neon_adjacency_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE))
dim(sub.data)
#Making sure that elements of dom.mat are read as numeric
mode(dom.mat) <- "numeric"

######################################
# Spatial Subset
######################################
# #Further subsetting while Balancing the unique number of species 
# Let's lower the resolution of the cells
#resolution: https://cran.r-project.org/web/packages/dggridR/vignettes/dggridR.html 
r = 8 #lower to include less observations // increase to include more observations
#generate hexagonal grid
dggs <- dgconstruct(res = r)
#> Resolution: 13, Area (km^2): 31.9926151554038, Spacing (km): 5.58632116604266, CLS (km): 6.38233997895802
#read in data
mod.data.cell <- sub.data %>%
  #get hexagonal cell id and week number for each checklist
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)
#sample from each instance
set.seed(20260226)
sub.data.sub1 <- mod.data.cell %>%
  group_by(species, domain_id, cell) %>% #unique species-domain instance -- balances species & domain, actual observations may vary
  sample_n(size = 1) %>%
  ungroup()
#subset data dimensions
dim(sub.data.sub1)

# Manual checks to ensure data look correct after subsetting 
a = sub.data.sub1 %>% count(species)
print(a, n=62)
b = sub.data.sub1 %>% count(species, domain_id)
print(b, n=677)

sub.data = sub.data.sub1

######################################
# Model calibration/validation groups
######################################
message("1) Create calibration and validation groups")
#Dummy data for testing code [remove slicing for actual runs]
# sub.data <- sub.data %>%
#   slice_sample(n = 1000)   # keep exactly n rows
#Training & testing groups -- for model validation
set.seed(run)   # Set a seed value for reproducibility purposes
#Randomly select 70% of the data (rows)
index <- sample(1:nrow(sub.data), size=floor(.70*nrow(sub.data)))
#Training group (initial model)
train.dat <- sub.data[index, ]
#Testing group
test.dat  <- sub.data[-index, ]
dim(train.dat); dim(test.dat)
#save space
rm(sub.data)
#Training data without woody species (for second model)
train.dat1 = train.dat[train.dat$funct_type == "non-woody", ]
test.dat1 = test.dat[test.dat$funct_type == "non-woody", ]
# Save training and test data for subsequent plotting 
save(train.dat1, 'trainingData1.Rdata')
save(test.dat1, 'testingData1.Rdata')

#####################
# Linear assumptions for Latitudinal Constraints Hypothesis Model (Functional type, but not nativity)
#####################
message("Exporting linear assumptions")

#Model with functional type
pdf(file = paste0("brms_fit_LinearAssumptions_Run", run, "_functionalType.pdf"))
mod <- lm(pdiff ~ rdiff*latitude + dispersal + func_type, data = train.dat1)
#1. Linearity & Homoscedasticity
par(mfrow = c(2,2))
plot(mod)
#2. Normality of residuals
hist(resid(mod), breaks = 100, main = "Histogram of Residuals", xlab = "Residuals")
#3. Q-Q plot
qqnorm(resid(mod)); qqline(resid(mod), col = "red")
dev.off()
#save space
rm(mod); gc()

x.time = system.time({ #START timer
  ###############################################################################
  #----- All woody & non-woody (w/o nativity status) :: With uncertainty --------
  ###############################################################################
  #Not including nativity status since only non-woody functional types are invasive
  message("4) Begin second brms model -- All woody & non-woody, w/o nativity")
  #Formua
  form <- bf(
    pdiff | se(psd, sigma = TRUE) ~
      me(rdiff, rsd)*latitude + funct_type + dispersal +
      (1 | species) + #random effect
      car(dom.mat, gr = domain_id)) #adjacency matrix
  message(form)
  #Setting priors
  priors <- c(
    set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
    set_prior("normal(0, 1)",         class = "b"),      # includes slope for me()
    set_prior("exponential(1)",       class = "sd"),     # RE SD
    set_prior("student_t(3, 0, 2.5)", class = "sigma")   # residual SD (beyond phen_sd)
  )
  #Fitting model
  fit <- brm(
    formula = form,
    data    = train.dat,
    data2 = list(dom.mat=dom.mat), # adjacency matrix data
    family  = gaussian(),           # identity link on [-1, 1]
    prior   = priors,
    chains  = 4, iter = 3000, warmup = 1000, cores = 40,
    threads = threading(160),
    seed    = run,
    control = list(adapt_delta = 0.95, max_treedepth = 13))
  message("Modeling complete! ... Saving model fit")
  saveRDS(fit, file = paste0("brms_fit_training_Run", run, "_funct_type.rds"))
  message("Predicting using the testing dataset")
  #predict testing data using the training model
  fit.t <- predict(fit, newdata = test.dat)
  saveRDS(fit.t, file = paste0("brms_fit_testing_Run", run, "_funct_type.rds"))
  
  #Save output
  message(cat("\tsaving output"))
  #Posterior Mean predictions
  #Posterior expected value per observation
  # train_fitted <- fitted(fit)  # matrix with Estimate, Est.Error, Q2.5, Q97.5
  # Posterior mean predictions
  train_fitted <- posterior_predict(fit, ndraws = 200) #subsampling or else it will crash R
  # train_mean <- train_fitted[, "Estimate"]
  train_mean <- colMeans(train_fitted)
  test_mean <- fit.t[, "Estimate"]  # if predict() returns a matrix with Estimate
  rmse_train <- sqrt(mean((train.dat$pdiff - train_mean)^2)); rmse_test <- sqrt(mean((test.dat$pdiff - test_mean)^2))
  r2.train <- cor(train.dat$pdiff, train_mean)^2; r2.test  <- cor(test.dat$pdiff, test_mean)^2
  # r2.train = bayes_R2(fit); r2.test = bayes_R2(fit.t)
  
  sink(paste0("brms_fit_Run", run, "_funct_type_output.txt"))
  cat("\nRMSE\n")
  cat("training: ", round(rmse_train, 4)); cat("\ntesting: ", round(rmse_test, 4))
  cat("\nR^2\n")
  cat("\ntraining: ", r2.train); cat("\ntesting: ", r2.test)
  cat("\n")
  cat("\nTraining\n")
  cat("Summary:\n"); print(summary(fit))
  cat("\nPrior summary:\n"); print(prior_summary(fit))
  cat("\nLook for the coefficient named b_me... — that's the slope for the latent RANGE.\n")
  print(fit)
  cat("\nTesting\n")
  cat("Summary:\n"); print(summary(fit.t))
  print(fit.t)
  sink()
  #-------------------
}) #END timer

sink(file = paste0("brms_fit_Run", run, "_funct_type_TotalRunTime.txt"))
min = round(unname(x.time[3]/60))
if(min > 59){
  hr = round(min/60)
  print(paste("Total run time: ", hr, " hours"))
  if(hr > 23){
    day = round(hr/24)
    print(paste("Total run time: ", day, " days"))
  }
} else{
  print(paste("Total run time: ", min, " minutes"))
}
sink()

#####################
# Linear assumptions of nativity model (Species Origina Hypothesis)
#####################
message("Exporting linear assumptions")

#Model with nativity (No functional type)
pdf(file = paste0("brms_fit_LinearAssumptions_Run", run, "_nativity.pdf"))
mod <- lm(pdiff ~ rdiff*latitude + rdiff*nativity + dispersal, data = train.dat1)
#1. Linearity & Homoscedasticity
par(mfrow = c(2,2))
plot(mod)
#2. Normality of residuals
hist(resid(mod), breaks = 100, main = "Histogram of Residuals", xlab = "Residuals")
#3. Q-Q plot
qqnorm(resid(mod)); qqline(resid(mod), col = "red")
dev.off()
#save space
rm(mod); gc()



###################
## Fit BRMS model:
###################
##   - Gaussian identity
##   - Response measurement error via  | se(phen_sd (pdiff), sigma = TRUE)
##   - Predictor measurement error via me(range_mean (rdiff), range_sd (rsd))
##   - Random intercept by species (can easily change to site or other factors)



x.time = system.time({ #START timer
  ###############################################################################
  #----------- All non-woody (w/ nativity status) :: With uncertainty -----------
  ###############################################################################
  #Removing woody species (and therefor function type) from the model (only non-woody invasive spp.)
  message("3) Begin first brms model -- All non-woody, w/ nativity; w/ uncertainty")
  #Formula
  form <- bf(
    pdiff | se(psd, sigma = TRUE) ~
      me(rdiff, rsd)*latitude + me(rdiff, rsd)*nativity + dispersal +
      (1 | species) + #random effect
      car(dom.mat, gr = domain_id)) #adjacency matrix
  message(form)
  #Setting priors
  priors <- c(
    set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
    set_prior("normal(0, 1)",         class = "b"),      # includes slope for me()
    set_prior("exponential(1)",       class = "sd"),     # RE SD
    set_prior("student_t(3, 0, 2.5)", class = "sigma")   # residual SD (beyond phen_sd)
  )
  #Fit model
  fit <- brm(
    formula = form,
    data    = train.dat1,
    data2 = list(dom.mat=dom.mat), # adjacency matrix data
    family  = gaussian(),           # identity link on [-1, 1]
    prior   = priors,
    chains  = 4, iter = 3000, warmup = 1000, cores = 40,
    threads = threading(160),
    seed    = run,
    control = list(adapt_delta = 0.95, max_treedepth = 13)
  )
  message(cat("\tModeling complete! ... Saving model fit"))
  saveRDS(fit, file = paste0("brms_fit_training_Run", run, "_nativity.rds"))
  message(cat("\tPredicting using the testing dataset"))
  #predict testing data using the training model
  fit.t <- predict(fit, newdata = test.dat1)
  saveRDS(fit.t, file = paste0("brms_fit_testing_Run", run, "_nativity.rds"))
  
  
  #Save output
  message(cat("\tsaving output"))
  #Posterior Mean predictions
  #Posterior expected value per observation
  # train_fitted <- fitted(fit)  # matrix with Estimate, Est.Error, Q2.5, Q97.5
  # Posterior mean predictions
  train_fitted <- posterior_predict(fit, ndraws = 200) #subsampling or else it will crash R
  # train_mean <- train_fitted[, "Estimate"]
  train_mean <- colMeans(train_fitted)
  test_mean <- fit.t[, "Estimate"]  # if predict() returns a matrix with Estimate
  rmse_train <- sqrt(mean((train.dat1$pdiff - train_mean)^2)); rmse_test <- sqrt(mean((test.dat1$pdiff - test_mean)^2))
  r2.train <- cor(train.dat1$pdiff, train_mean)^2; r2.test  <- cor(test.dat1$pdiff, test_mean)^2
  # r2.train = bayes_R2(fit); r2.test = bayes_R2(fit.t)
  
  sink(paste0("brms_fit_Run", run, "_nativity_output.txt"))
  cat("\nRMSE\n")
  cat("training: ", round(rmse_train, 4)); cat("\ntesting: ", round(rmse_test, 4))
  cat("\nR^2\n")
  cat("\ntraining: ", r2.train); cat("\ntesting: ", r2.test)
  cat("\n")
  cat("\nTraining\n")
  cat("Summary:\n"); print(summary(fit))
  cat("\nPrior summary:\n"); print(prior_summary(fit))
  cat("\nLook for the coefficient named b_me... — that's the slope for the latent RANGE.\n")
  print(fit)
  cat("\nTesting\n")
  cat("Summary:\n"); print(summary(fit.t))
  print(fit.t)
  sink()
  #-------------------
}) #END timer

sink(file = paste0("brms_fit_Run", run, "_nativity_TotalRunTime.txt"))
min = round(unname(x.time[3]/60))
if(min > 59){
  hr = round(min/60)
  print(paste("Total run time: ", hr, " hours"))
  if(hr > 23){
    day = round(hr/24)
    print(paste("Total run time: ", day, " days"))
  }
} else{
  print(paste("Total run time: ", min, " minutes"))
}
sink()






# #----------- All non-woody (w/ nativity status) :: No priors -----------
# message("4) Begin second brms model -- All non-woody, w/ nativity; NO PRIORS")
# #Removing woody species (and therefor function type) from the model (only non-woody invasive spp.)
# #Formula
# form <- bf(
#   pdiff  ~
#     rdiff + latitude + nativity + dispersal + rdiff*latitude +
#     (1 | species) + #random effect
#     car(dom.mat, gr = domain_id)) #adjacency matrix
# message(form)
# 
# #Fitt model
# fit <- brm(
#   formula = form,
#   data    = train.dat1,
#   data2 = list(dom.mat=dom.mat), # adjacency matrix data
#   family  = gaussian(),           # identity link on [-1, 1]
#   #prior   = priors,
#   chains  = 4, iter = 3000, warmup = 1000, cores = 4,
#   #threads = threading(16),
#   seed    = 20250909,
#   #control = list(adapt_delta = 0.95, max_treedepth = 13)
# )
# message(cat("\tModeling complete! ... Saving model fit"))
# saveRDS(fit, file = file.path(L2, paste0("brms_fit_Run", run, "_training_nativity_NP.rds")))
# message(cat("\tPredicting using the testing dataset"))
# #predict testing data using the training model
# fit.t <- predict(fit, newdata = test.dat1)
# saveRDS(fit.t, file = file.path(L2, paste0("brms_fit_Run", run, "_testing_nativity_NP.rds")))
# 
# }) #END timer
# 
# sink(file = file.path(graphs, paste0("brms_fit_Run", run, "_nativity_TotalRunTime.txt")))
# min = round(x.time[3]/60)
# if(min > 59){
#   hr = round(min/60)
#   message("Total run time: ", hr, " hours")
#   if(hr > 23){
#     day = round(hr/24)
#     message("Total run time: ", day, " days")
#   }
# } else{
#   message("Total run time: ", min, " minutes")
# }
# sink()
# 
# 
# #--- save output ---
# message(cat("\tsaving output"))
# #Posterior Mean predictions
# #Posterior expected value per observation
# # train_fitted <- fitted(fit)  # matrix with Estimate, Est.Error, Q2.5, Q97.5
# # Posterior mean predictions
# train_fitted <- posterior_predict(fit, ndraws = 200)
# # train_mean <- train_fitted[, "Estimate"]
# train_mean <- colMeans(train_fitted)
# test_mean <- fit.t[, "Estimate"]  # if predict() returns a matrix with Estimate
# rmse_train <- sqrt(mean((train.dat1$pdiff - train_mean)^2)); rmse_test <- sqrt(mean((test.dat1$pdiff - test_mean)^2))
# r2.train <- cor(train.dat1$pdiff, train_mean)^2; r2.test  <- cor(test.dat1$pdiff, test_mean)^2
# # r2.train = bayes_R2(fit); r2.test = bayes_R2(fit.t)
# 
# sink(file.path(graphs, paste0("brms_fit_Run", run, "_nativity_output_NP.txt")))
# cat("\nRMSE\n")
# cat("training: ", round(rmse_train, 4)); cat("\ntesting: ", round(rmse_test, 4))
# cat("\nR^2\n")
# cat("\ntraining: ", r2.train); cat("\ntesting: ", r2.test)
# cat("\n")
# cat("\nTraining\n")
# cat("Summary:\n"); print(summary(fit))
# cat("\nPrior summary:\n"); print(prior_summary(fit))
# cat("\nLook for the coefficient named b_me... — that's the slope for the latent RANGE.\n")
# print(fit)
# cat("\nTesting\n")
# cat("Summary:\n"); print(summary(fit.t))
# print(fit.t)
# sink()
# #-------------------