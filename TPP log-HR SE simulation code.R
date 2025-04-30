# © Copyright reserved by the Authors and PharmaQuant Insights Pvt. Ltd, India @2025
# Created by Abhirup Dutta Majumdar
# Quality checked by Sekhar Kumar Dutta 
# Note: The user should always fully credit the authors of this paper in every private and public work conducted. If ever found in violation, there could be legal and copyright infringement consequences
##### Estimation of standard error######
setwd("") #add your directory
library(tidyverse)
library(survival)
library(survminer)
library(flexsurv)
library(simsurv)
library(ggplot2)
library(progress)
library(MASS)
rm(list = ls())

#Loading reconstructed km data from Guyot 2012 output
sunitnib_pfs<-read.csv("IPD_CheckMate 214_PFS_Sunitinib.csv", header = T) # add your path to the reconstructed file for anchor study comparator arm

########## RP Model fitting ###########
### One internal knot ###
flexph1<-flexsurvspline(Surv(time, death) ~ 1, 
                        data = sunitnib_pfs, k = 1)
BIC(flexph1)
### Two internal knots ### 
flexph2 <- flexsurvspline(Surv(time, death) ~ 1, 
                          data = sunitnib_pfs, k = 2)
flexph2
BIC(flexph2)
### Three internal knots ### 
flexph3<- flexsurvspline(Surv(time, death) ~ 1, 
                         data = sunitnib_pfs, k = 3)
flexph3
BIC(flexph3)
### Four internal knots ### 
flexph4<-flexsurvspline(Surv(time, death) ~ 1, 
                        data = sunitnib_pfs, k = 4)
flexph4
BIC (flexph4)

### Five internal knots ### 
flexph5<-flexsurvspline(Surv(time, death) ~ 1, 
                        data = sunitnib_pfs, k = 5)
flexph5
BIC (flexph5)

### Six internal knots ### 
flexph6<-flexsurvspline(Surv(time, death) ~ 1, 
                        data = sunitnib_pfs, k = 6)
flexph6
BIC (flexph6)

##########Plotting of fitted curves over KM curve###########
par(mfrow = c(2, 3), cex = 0.55) # graphics parameters
plot(flexph1, 
     main = "Flexible par model knot1 ", 
     ylab = "Survival probability", 
     xlab = "Time")
plot(flexph2, 
     main = "Flexible par model knot2", 
     ylab = "Survival probability", 
     xlab = "Time")
plot(flexph3, 
     main = "Flexible par model knot3", 
     ylab = "Survival probability", 
     xlab = "Time")
plot(flexph4, 
     main = "Flexible par model knot4", 
     ylab = "Survival probability", 
     xlab = "Time")
plot(flexph5, 
     main = "Flexible par model knot5", 
     ylab = "Survival probability", 
     xlab = "Time")
plot(flexph6, 
     main = "Flexible par model knot6", 
     ylab = "Survival probability", 
     xlab = "Time")


########## Matrix of BIC of all fitted models#############
BIC<- matrix(c(BIC(flexph1),BIC(flexph2), BIC(flexph3), BIC(flexph4),BIC(flexph5), BIC(flexph6),
             nrow = 6,ncol = 1, byrow = T)
rownames(BIC)<- c("Flex CS knot1",
                      "Flex CS knot2","Flex CS knot3","Flex CS knot4","Flex CS knot5","Flex CS knot6" )
colnames(BIC)<- c("BIC")
 
BIC_min<- rownames(BIC)[which.min(BIC[, 1])] ##CS2
print(BIC_min) 

######Simulating parametric standard errors based on best-fit#####

set.seed(123)  # Set a single master seed for the bootstrap as per Morris 2019 (the users may wish to save the RNG state if there are convergence issues, not showed here)

# Define custom HR values
HR_values <- c(0.1, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25, 0.275, 0.30, 0.325, 0.35, 0.375, 0.40,0.425, 0.45, 0.475, 0.50, 0.525, 0.55, 0.575, 0.60, 0.625, 0.65, 
               0.675, 0.70, 0.725, 0.75, 0.775, 0.80, 0.825, 0.85, 0.875, 0.90, 0.925, 0.95, 0.975, 1.00)
# Save HR values to a CSV file
write.csv(data.frame(HR_values), "HR_values.csv", row.names = FALSE)
cat("\n HR_values saved to HR_values.csv\n")

# Number of simulations (targeting 95% coverage and a Monte Carlo SE of bias value to be lower than 0.005, as per Morris 2019)
nsim <- 1600 

# Start timer
start_time <- Sys.time()

# Function to update anch_mod for each HR
update_model <- function(HR) {
  flexph2 <- flexph2  
  
  # Create betas based on HR
  betas1 <- data.frame(
    gamma0 = rep(flexph2$res.t["gamma0", "est"], 1096), 
    gamma1 = rep(flexph2$res.t["gamma1", "est"], 1096), 
    gamma2 = rep(flexph2$res.t["gamma2", "est"], 1096), 
    gamma3 = rep(flexph2$res.t["gamma3", "est"], 1096),
    arm = rep(log(HR), 1096)  
  )
  
  return(list(flexph2 = flexph2, betas1 = betas1))
}

# Log cumulative hazard function
logcumhaz <- function(t, x, betas, knots) {
  basis <- flexsurv::basis(knots, log(t))
  res <- 
    betas[["gamma0"]] * basis[[1]] + 
    betas[["gamma1"]] * basis[[2]] +
    betas[["gamma2"]] * basis[[3]] +
    betas[["gamma3"]] * basis[[4]] +
    betas[["arm"]] * x[["arm"]] # modelled assuming proportional-hazards under the RP structure
  res
}

# Simulation function
sim_run <- function(true_mod, betas1, HR) {
  # Covariate data
  cov <- data.frame(id = 1:1096, arm = rbinom(1096, 1, 0.5))
  
  # Simulate event times
  dat <- simsurv(betas = betas1, 
                 x = cov, 
                 knots = true_mod$flexph2$knots, 
                 logcumhazard = logcumhaz, 
                 maxt = 66, # max time of follow-up from Checkmate-214 
                 interval = c(1E-8, 100000))
  
  # Merge covariate data and event times
  dat <- merge(cov, dat)
  
  # Fit flexible survival model
  flex_mod <- flexsurvspline(Surv(eventtime, status) ~ arm, data = dat, k = 2) #need to change based on best-fit
  
  # Compute log HR and confidence interval
  true_loghr <- log(HR)
  flex_loghr <- flex_mod$coefficients[["arm"]]
  ses <- flex_mod$res.t["arm", "se"]
  cil <- flex_loghr + qnorm(.025) * ses
  ciu <- flex_loghr + qnorm(.975) * ses
  
  # Return bias and coverage results
  return(c(true_loghr = true_loghr, 
           flex_bias = flex_loghr - true_loghr, 
           se_flex = ses, 
           events = flex_mod$events, 
           coverage = (true_loghr > cil) && (true_loghr < ciu)))
}

# Run simulations for each HR value
all_results <- list()

for (HR in HR_values) {
  # Update the model for the current HR
  model_data <- update_model(HR)
  
  # Progress bar
  pb <- progress_bar$new(
    format = paste0("HR = ", HR, " [:bar] :percent in :elapsed"),
    total = nsim,
    width = 60
  )
  
  # Run simulations
  results <- vector("list", nsim)
  for (i in 1:nsim) {
    results[[i]] <- sim_run(true_mod = model_data, betas1 = model_data$betas1, HR = HR)
    pb$tick()
  }
  
  # Convert results to DataFrame and add HR column
  results_df <- as.data.frame(do.call(rbind, results))
  results_df$HR <- HR  # Add HR column
  
  # Save to CSV file
  filename <- sprintf("simulation_results_HR_%.2f.csv", HR)
  write.csv(results_df, filename, row.names = FALSE)
  cat(sprintf("\n saved: %s\n", filename))  # Print confirmation
  
  # Store results in list for optional combined file
  all_results[[as.character(HR)]] <- results_df
}

# Combine all results into a single data frame (optional)
final_results <- do.call(rbind, all_results)

# Save a single file with all results (optional)
write.csv(final_results, "all_simulation_results_pfS.csv", row.names = FALSE)

# End timer
end_time <- Sys.time()
print(end_time - start_time)

cat("\n All simulations completed successfully \n")
