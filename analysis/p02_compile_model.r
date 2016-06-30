sink("./analysis/output/l02_compile_model.rlog", split=TRUE)
## p03_fit_models_cluster.r
## Fit Poisson models for each cancer (using the Imperial cluster)


## required packages
req <- c("rstan")
lapply(req, library, character.only = TRUE)
rm(req)

set.seed(893749)
sessionInfo()
##############################################

poisson_ods_full_centred <- stan_model("./analysis/m01_poisson_ods_full_centred.stan",
                                       verbose=TRUE)

save(poisson_ods_full_centred,
     file="./analysis/m01_poisson_ods_full_centred_compiled.Rdta")

##############################################
sink()
