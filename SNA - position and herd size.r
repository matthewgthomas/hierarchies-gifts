##
## does position in gift network predict herd size and/or no. kids?
##
library(tidyverse)
library(rstanarm)
library(loo)
library(bayesplot)
library(futile.logger)

options(mc.cores = parallel::detectCores())

source("init if necessary.r")
source("model comparison functions.r")

# data for analyses
d.cattle = people %>% 
  select(HerderID, Pop, cattle.z, Age.z, Sexb, gift.deg.in, gift.betweenness, gift.eigen) %>% 
  na.omit(.)

# uncomment this if you want to fit the models to only the Saami sites
# d.cattle = d.cattle %>% 
#   filter(Pop %in% c("karasjok", "kautokeino")) %>% 
#   mutate(Pop = forcats::fct_drop(Pop))


################################################################
## Plot descriptives
##
# plot association between age and herd size (with LOESS smoothers)
people %>% 
  select(Pop, Age, HerdSize) %>% 
  na.omit(.) %>% 
  mutate(PopName = convert_pop(Pop)) %>% 
  
  ggplot(aes(x=Age, y=HerdSize)) + 
    geom_point(aes(colour=PopName), alpha=0.3) + 
    geom_smooth(aes(colour=PopName, fill=PopName), se=T) +
    facet_wrap(~PopName, scales="free_y") +
    ylab("Herd size") +
    common_theme

ggsave(filename=file.path(plots.dir, "Age and herd size.png"), width=250, height=150, units="mm")
ggsave(filename=file.path(plots.dir, "Age and herd size.pdf"), width=200, height=120, units="mm")

# does age predict herd size? No.
# library(broom)
# 
# d.cattle %>% 
#   group_by(Pop) %>% 
#   #do(tidy(lm(HerdSize ~ Age, data=.))) %>% 
#   do(tidy(cor.test(~ HerdSize + Age, data=.))) %>% 
#   mutate(p.value = round(p.value, 4))


################################################################
## Predicting herd size
##
flog.info("Fitting herd size models")

m.herd.null = stan_lmer(cattle.z ~ 1 + (1 | Pop),
                        data=d.cattle,
                        prior_intercept=normal(0, 5), 
                        adapt_delta=0.99, chains=4, cores=2)

m.herd.control = stan_lmer(cattle.z ~ Age.z + Sexb + (1 | Pop),
                           data=d.cattle,
                           prior_intercept=normal(0, 5), prior=normal(0, 5), 
                           adapt_delta=0.99, chains=4, cores=2)

m.herd.deg = stan_lmer(cattle.z ~ gift.deg.in + Age.z + Sexb + (1 | Pop),
                       data=d.cattle,
                       prior_intercept=normal(0, 5), prior=normal(0, 5),
                       adapt_delta=0.99, chains=4, cores=2)

m.herd.deg_int = stan_lmer(cattle.z ~ gift.deg.in * Age.z + Sexb + (1 | Pop),
                           data=d.cattle,
                           prior_intercept=normal(0, 5), prior=normal(0, 5),
                           adapt_delta=0.99, chains=4, cores=2)

m.herd.btwn = stan_lmer(cattle.z ~ gift.betweenness + Age.z + Sexb + (1 | Pop),
                        data=d.cattle,
                        prior_intercept=normal(0, 5), prior=normal(0, 5),
                        adapt_delta=0.99, chains=4, cores=2)

m.herd.eigen = stan_lmer(cattle.z ~ gift.eigen + Age.z + Sexb + (1 | Pop),
                         data=d.cattle,
                         prior_intercept=normal(0, 5), prior=normal(0, 5),
                         adapt_delta=0.99, chains=4, cores=2)


################################################################
## Model selection
##
flog.info("Comparing herd size models")

loo.herd.null    = loo(m.herd.null)
loo.herd.control = loo(m.herd.control)
loo.herd.deg     = loo(m.herd.deg)
loo.herd.btwn    = loo(m.herd.btwn)
loo.herd.eigen   = loo(m.herd.eigen)

loo.herd.deg_int = loo(m.herd.deg_int)
ll.herd.deg_int = log_lik(m.herd.deg_int)

ll.herd.null    = log_lik(m.herd.null)
ll.herd.control = log_lik(m.herd.control)
ll.herd.deg     = log_lik(m.herd.deg)
ll.herd.btwn    = log_lik(m.herd.btwn)
ll.herd.eigen   = log_lik(m.herd.eigen)

m.herd.loos         = compare_models_delta(loo.herd.null, loo.herd.control, loo.herd.deg, loo.herd.deg_int, loo.herd.btwn, loo.herd.eigen)
m.herd.loos.weights = model_weights(list(ll.herd.null, ll.herd.control, ll.herd.deg, ll.herd.deg_int, ll.herd.btwn, ll.herd.eigen))
m.herd.loos.sel     = model_select(list(ll.herd.null, ll.herd.control, ll.herd.deg, ll.herd.deg_int, ll.herd.btwn, ll.herd.eigen))

##
## save model comparison table
##
m_names = c(m.herd.null    ="(Varying intercepts only)", 
            m.herd.control = "Controls only: Age + Sex", 
            m.herd.deg     = "In-degree", 
            m.herd.deg_int = "In-degree × Age", 
            m.herd.btwn    = "Betweenness", 
            m.herd.eigen   = "Eigenvector centrality")

names(m.herd.loos.weights) = c("m.herd.null", "m.herd.control", "m.herd.deg", "m.herd.deg_int", "m.herd.btwn", "m.herd.eigen")

# add stacking weights to table (in correct order)
m.herd.loos$stacking = round(m.herd.loos.weights[row.names(m.herd.loos)], 4)  # m.herd.loos.weights[ m_order ]

# add model names to table
#m_order = as.integer( gsub("model([0-9])", "\\1", row.names(m.herd.loos)) )  # get order of models from lowest to highest delta
row.names(m.herd.loos) = m_names[row.names(m.herd.loos)]  #m_names[ m_order ]

# calculate and save variance partition coefficients
vpc.herd.deg = VarCorr(m.herd.deg) %>% 
  dplyr::as_data_frame() %>% 
  mutate(vpc = vcov / sum(vcov))

##
## save everything
##
write.csv( subset(m.herd.loos, select=c(delta, d_se, stacking)), 
           file.path(results.dir, "social network analysis of herd size.csv"), row.names=T)

write.csv(vpc.herd.deg, file.path(results.dir, "herd size model VPCs.csv"), row.names=F)

save(d.cattle, 
     m.herd.null, m.herd.control, m.herd.deg, m.herd.btwn, m.herd.eigen, 
     loo.herd.null, loo.herd.control, loo.herd.deg, loo.herd.deg_int, loo.herd.btwn, loo.herd.eigen,
     ll.herd.null, ll.herd.control, ll.herd.deg, ll.herd.deg_int, ll.herd.btwn, ll.herd.eigen,
     m.herd.loos, m.herd.loos.weights, m.herd.loos.sel, vpc.herd.deg,
     file=file.path(models.dir, "herd size models.RData"))



########################################
## plot parameter estimates
##
flog.info("Plotting parameter estimates")

draws = as.data.frame(m.herd.deg)
names(draws) = c("(Intercept)", "No. gifts received", "Age (z-score)", "Sex (0 = male)",
                 "Tibet: Cairima intercept", "Tibet: Doulong intercept", "Tibet: Jilehe intercept",
                 "Finnmark: Karasjok intercept", "Finnmark: Kautokeino intercept", "Tibet: Tawa intercept",
                 "Between-herder variance", "Between-site variance")

# put param names in a better order (mainly study sites)
draws = draws %>% 
  select("(Intercept)", "No. gifts received", "Age (z-score)", "Sex (0 = male)", 
         starts_with("Finnmark"), starts_with("Tibet"), everything())
  
mcmc_intervals(draws, prob = 0.8, prob_outer = 0.95) +
  xlab("Coefficient") +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines") 
        ,strip.background = element_blank()
        ,strip.placement = "outside"
        ,plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.background = element_blank()
  )

ggsave(file.path(plots.dir, "Prediction intervals - herd sizes from gifts received.png"), height=100, width=200, units="mm")
ggsave(file.path(plots.dir, "Prediction intervals - herd sizes from gifts received.pdf"), height=100, width=200, units="mm")


########################################
## predicted deviation in herd size given no. gifts received
##
flog.info("Plotting predictions")

new.data = expand.grid(
  gift.deg.in = seq(from=min(d.cattle$gift.deg.in), to=max(d.cattle$gift.deg.in)),
  Age.z = 0,
  Sexb = 0:1,
  Pop = "",
  cattle.z = 0
)

post = posterior_predict(m.herd.deg, new.data, re.form = NA)
pred = posterior_linpred(m.herd.deg, new.data, transform=T, re.form = NA)

quants = apply(post, 2, quantile, probs = c(0.025, 0.5, 0.975))  # quantiles over mcmc samples
quants2 = apply(pred, 2, quantile, probs = c(0.025, 0.5, 0.975))  # quantiles over mcmc samples
row.names(quants) = c("sim.lwr", "sim.med", "sim.upr")
row.names(quants2) = c("lwr", "herd.pred", "upr")

new.data = cbind(new.data, t(quants), t(quants2))

ggplot(new.data, aes(x=gift.deg.in, y=herd.pred)) +
  geom_line(aes(colour=factor(Sexb))) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=factor(Sexb)), alpha=0.2) +

  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  
  xlab("No. gifts received") +
  ylab("Herd size (z-score)") +
  
  # plot relatedness on log scale with breaks at familiar points
  scale_x_continuous(breaks=c(0,2,4,6,8,10)) +

  common_theme

ggsave(filename=file.path(plots.dir, "Predicted herd size from gifts received.png"), width=100, height=80, units="mm")
ggsave(filename=file.path(plots.dir, "Predicted herd size from gifts received.pdf"), width=100, height=80, units="mm")
