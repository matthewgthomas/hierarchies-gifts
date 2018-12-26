##
## calculate LOOs for each of the multilevel models
##
library(rstanarm)
library(loo)

source("init if necessary.r")
source("model comparison functions.r")  # load customised functions for model comparison

# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())

load(file.path(models.dir, gift_models_file))  # run 'multilevel models - fit.r' if this doesn't exist


#######################################################################
## Model comparison and diagnostics
##
flog.info("Calculating LOOs")

# calculate loos (warning: this takes many hours)
loo.null.nest   = loo_sparse(mlm.null.nested);  flog.info("... finished null")
loo.r.re        = loo_sparse(mlm.r.re)       ;  flog.info("... finished r")
loo.siida.re    = loo_sparse(mlm.siida.re)   ;  flog.info("... finished siida")
loo.r_siida.fe  = loo_sparse(mlm.r_siida.fe) ;  flog.info("... finished r + siida - intercepts")
loo.r_siida.re  = loo_sparse(mlm.r_siida.re) ;  flog.info("... finished r + siida - slopes")
loo.full.fe     = loo_sparse(mlm.full.fe)    ;  flog.info("... finished r x siida - intercepts")
loo.full        = loo_sparse(mlm.full)       ;  flog.info("... finished r x siida - slopes")

# set model names in the loo objects
attr(loo.null.nest , "name") = "mlm.null.nested"
attr(loo.r.re,       "name") = "mlm.r.re"
attr(loo.siida.re,   "name") = "mlm.siida.re"
attr(loo.r_siida.fe, "name") = "mlm.r_siida.fe"
attr(loo.r_siida.re, "name") = "mlm.r_siida.re"
attr(loo.full.fe,    "name") = "mlm.full.fe"
attr(loo.full,       "name") = "mlm.full"

# calculate deltas for models
# (don't trust these weights -- use the model_weights.loo() function, below)
loo_comp = compare_models_delta(loo.null.nest, loo.r.re, loo.siida.re, loo.r_siida.fe, loo.r_siida.re, loo.full.fe, loo.full)

loos = list(loo.null.nest, loo.r.re, loo.siida.re, loo.r_siida.fe, loo.r_siida.re, loo.full.fe, loo.full)

##
## model comparison using stacking and pseudo-Bayesian model averaging
##
flog.info("Calculating model weights")

mlm_weights        = model_weights.loo(loos, method="stacking")
# mlm_weights.bma    = model_weights.loo(loos, method="pseudobma", BB=F)  # this puts all the weight on model 7
# mlm_weights.bma_bb = model_weights.loo(loos, method="pseudobma", BB=T)  # so does this

# put model names into weights arrays
loos.names = unlist(lapply(loos, function(l) attr(l, "model")))
names(mlm_weights) = loos.names
# names(mlm_weights.bma) = loos.names
# names(mlm_weights.bma_bb) = loos.names

flog.info("Saving model comparisons")
save(loos, loo.null.nest, loo.r.re, loo.siida.re, loo.r_siida.fe, loo.r_siida.re, loo.full.fe, loo.full, loo_comp, mlm_weights,
     file=comparison_data_file)
