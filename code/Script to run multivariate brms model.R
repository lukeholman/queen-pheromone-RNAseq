library(brms)
library(purrr)
library(dplyr)
library(parallel)

# Note: make sure to have the 'eigen.data' object from the main .rmd script loaded into R

# Define the 5 models' formulae
model_formulae <- paste(map2(rep("cbind(m1,m2,m3,m4,m5,m6,m7,m8,m9) ~", 5),
                                    c("Species * Treatment", "Species + Treatment", "Treatment", "Species", "1"), ~ paste(.x, .y)), "+ (1 | p | colony)")

# Make the variance in the data equal to 1 for each eigengene (it's equal across modules already, but not equal to 1)
for(i in 1:9) eigen.data$Eigengene[eigen.data$Module == paste("Module", i)] <- as.numeric(scale(eigen.data$Eigengene[eigen.data$Module == paste("Module", i)]))

# Define a function to run the same model on each formula
run_eigen_model <- function(formula){
  if(!grepl("~ 1", formula)){ # If NOT doing the intercept-only model, fit weak priors on the fixed effects
    brm(as.formula(formula), 
        data = eigen.data %>% mutate(Module = gsub("Module ", "m", Module)) %>% select(-id) %>% spread(Module, Eigengene), 
        family = "gaussian",
        seed = 1,   # reproducible results
        cores = 1, chains = 4, iter = 5000, # Run 4 chains, in parallel, for 8000 iterations (3000 of which are discarded as burn-in)
        control = list(adapt_delta = 0.9999, max_treedepth = 15), # Slower results, but more reliable
        prior = set_prior("normal(0,5)", class = "b"), # Specify a weakly informative, normal prior for all fixed effects (units are Cohen's d)
        save_all_pars = TRUE) # needed to calculate model posterior probabilities
  }
  else {
    brm(as.formula(formula),  # No priors are needed for the intercept-only model
        data = eigen.data %>% mutate(Module = gsub("Module ", "m", Module)) %>% select(-id) %>% spread(Module, Eigengene), 
        family = "gaussian",
        seed = 1,   # reproducible results
        cores = 1, chains = 4, iter = 5000, # Run 4 chains, in parallel, for 8000 iterations (3000 of which are discarded as burn-in)
        control = list(adapt_delta = 0.9999, max_treedepth = 15), # Slower results, but more reliable
        save_all_pars = TRUE) # needed to calculate model posterior probabilities
  }

}

# Run 5 models on the 5 formula, in parallel using 5 cores
model_list <- mclapply(model_formulae, run_eigen_model, mc.cores = 5)

post_model_probs <- round(post_prob(model_list[[1]], model_list[[2]], model_list[[3]], model_list[[4]], model_list[[5]]), 2)
names(post_model_probs) <- c("Treatment x Species", "Treatment + Species", "Treatment", "Species", "Intercept only")
post_model_probs <- melt(post_model_probs) %>% rename(`Posterior model probability` = value) 

# Save the best model summary, and the comparison table
saveRDS(summary(model_list[[3]]), file = "data/brms_model_summary.rds")
saveRDS(post_model_probs, file = "data/brms_model_comparisons.rds")

# Annotations for figure 4 inset
treatment_effects <- readRDS("data/brms_model_summary.rds")$fixed[10:18, ] %>%
  round(2)
treatment_effects <- unname(apply(treatment_effects, 1, function(x) paste("d = ", format(x[1], nSmall = 2), ", (", format(x[3], nSmall = 2), " to ", format(x[4], nSmall = 2), ")", sep = "")))
