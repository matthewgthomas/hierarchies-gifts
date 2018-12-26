##
## Fit Bayesian hierarchical models using rstanarm
##
## models including varying intercepts for study site (Pop) and egos nested in sites (Pop:Ego)
## - some models also include varying slopes for sites (but not egos - we're only interested in between-site variation, not between-person)
##
library(rstanarm)
options(mc.cores = parallel::detectCores())

source("init if necessary.r")

##
## intercept-only - nested
## 24246.9 seconds (Total)
##
mlm.null.nested = stan_glmer(GiftGiven ~ 1 + (1 | Pop/Ego.ID),
                             data=d.long, family=binomial,
                             prior_intercept=normal(0, 5),
                             chains=1, cores=1)
save(mlm.null.nested, file=file.path(models.dir, "mlm - varying intercepts only.RData"))

##
## relatedness only - varying slopes
## 86973.1 seconds (Total)
##
mlm.r.re = stan_glmer(GiftGiven ~ r + (r | Pop) + (1 | Pop:Ego.ID), 
                      data=d.long, family=binomial,
                      prior_intercept=normal(0, 5), prior=normal(0,2), prior_covariance=decov(regularization=2),
                      chains=1, cores=1)
save(mlm.r.re, file=file.path(models.dir, "mlm - r - varying intercepts and slopes.RData"))

##
## herding group only - varying slopes
## 47142.8 seconds (Total)
##
mlm.siida.re = stan_glmer(GiftGiven ~ SameHerdingGroup + (SameHerdingGroup | Pop) + (1 | Pop:Ego.ID), 
                          data=d.long, family=binomial,
                          prior_intercept=normal(0, 5), prior=normal(0,2), prior_covariance=decov(regularization=2),
                          chains=1, cores=1)
save(mlm.siida.re, file=file.path(models.dir, "mlm - siida - varying intercepts and slopes.RData"))

##
## r + herding group - varying intercepts
## 30491.6 seconds (Total)
##
mlm.r_siida.fe = stan_glmer(GiftGiven ~ r + SameHerdingGroup + (1 | Pop/Ego.ID), 
                            data=d.long, family=binomial,
                            prior_intercept=normal(0, 5), prior=normal(0,2), prior_covariance=decov(regularization=2),
                            QR=T, chains=1, cores=1)
save(mlm.r_siida.fe, file=file.path(models.dir, "mlm - r plus siida - varying intercepts.RData"))

##
## r + herding group - varying intercepts and slopes
## 112458 seconds (Total)
##
mlm.r_siida.re = stan_glmer(GiftGiven ~ r + SameHerdingGroup + (r + SameHerdingGroup | Pop) + (1 | Pop:Ego.ID), 
                            data=d.long, family=binomial,
                            prior_intercept=normal(0, 5), prior=normal(0,2), prior_covariance=decov(regularization=2),
                            QR=T, chains=1, cores=1)
save(mlm.r_siida.re, file=file.path(models.dir, "mlm - r plus siida - varying intercepts and slopes.RData"))

##
## r x herding group - varying intercepts by egos in populations
## 42053.4 seconds (Total)
##
mlm.full.fe = stan_glmer(GiftGiven ~ r * SameHerdingGroup + (1 | Pop/Ego.ID), 
                         data=d.long, family=binomial,
                         prior_intercept=normal(0, 5), prior=normal(0,2), prior_covariance=decov(regularization=2),
                         QR=T, chains=1, cores=1)
save(mlm.full.fe, file=file.path(models.dir, "mlm - r x siida - varying intercepts.RData"))

##
## r x herding group - varying intercepts and slopes by egos in populations
## 264341 seconds (Total)
##
mlm.full = stan_glmer(GiftGiven ~ r * SameHerdingGroup + (r * SameHerdingGroup | Pop) + (1 | Pop:Ego.ID), 
                      data=d.long, family=binomial,
                      prior_intercept=normal(0, 5), prior=normal(0,2), prior_covariance=decov(regularization=2),
                      QR=T, chains=1, cores=1)
save(mlm.full, file=file.path(models.dir, "mlm - r x siida - varying intercepts and slopes.RData"))

# save everything
save(d.long, mlm.null.nested, mlm.r.re, mlm.siida.re, mlm.r_siida.fe, mlm.r_siida.re, mlm.full.fe, mlm.full,
     file=gift_models_file)
