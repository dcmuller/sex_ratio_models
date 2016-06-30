#sink("./analysis/output/l03_fit_models_cluster.rlog", split=TRUE)
# can't use sink as will be called in parallel
## p03_fit_models_cluster.r
## Fit Poisson models for each cancer (using the Imperial cluster)

## get command line argument
args <- commandArgs(trailingOnly=TRUE)
idx <- as.numeric(args[1])
print(idx)

## required packages
req <- c("rstan")
lapply(req, library, character.only = TRUE)
rm(req)

set.seed(3749)
sessionInfo()
##############################################
load("./data/d00_ci5.Rdta")
load("./analysis/m01_poisson_ods_full_centred_compiled.Rdta")

## prepare data to pass to stan model
d <- ci5_bycancer[[idx]]
dlist <- with(d, list(N=nrow(d),
                      n_coef=1,
                      y=n_ca,
                      x=as.matrix(male),
                      reg=as.numeric(registry_cde),
                      n_reg=max(as.numeric(registry_cde)),
                      age=as.numeric(age_group),
                      n_age=max(as.numeric(age_group)),
                      year=as.numeric(factor(year)),
                      n_year=max(as.numeric(factor(year))),
                      log_offset=log(n_pop)))

fit <- sampling(poisson_ods_full_centred,
                data=dlist,
                iter=2000,
                chains=10, cores=10,
                control=list(max_treedepth=15))

fn <- paste0("./analysis/output/o03_fit_", gsub(" ", "_", names(ci5_bycancer)[idx]), ".Rdta")
save(fit, file=fn)
