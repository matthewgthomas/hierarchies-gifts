##
## Analyse varying intercepts/slopes models of gift giving
##
## run 'multilevel models - fit.r' first if you need to fit the models (warning: will take a LOOOOOooooOOOONG time)
## and run 'multilevel models - compare.r' to calculate information criteria etc. (takes about a day)
##
## install development version of loo for model_weights()
## > devtools::install_github("stan-dev/loo")
##
library(rstanarm)
library(loo)
library(bayesplot)
library(tidyverse)
library(xlsx)
library(RColorBrewer)

# load models and comparisons
if (file.exists(gift_models_file)) {
  load(gift_models_file)  # run `multilevel models - fit.r` if this doesn't exist
} else {
  stop("You need to fit the multilevel models before analysing them")
}

if (file.exists(comparison_data_file)) {
  load(comparison_data_file)  # run `multilevel models - compare.r` if this doesn't exist
} else {
  stop("You need to compare the multilevel models before analysing them")
}

# some useful lists
mlms = list(mlm.null.nested, mlm.r.re, mlm.siida.re, mlm.r_siida.fe, mlm.r_siida.re, mlm.full.fe, mlm.full)
models_list = c("mlm.null.nested" = "Null - nested", 
                "mlm.r.re" = "r", 
                "mlm.siida.re" = "siida",
                "mlm.r_siida.fe" = "r + siida", 
                "mlm.r_siida.re" = "r + siida - slopes", 
                "mlm.full.fe" = "r * siida", 
                "mlm.full" = "r * siida - slopes")


#######################################################################
## Check model selection table
##
loo_comp
round(mlm_weights, 2)
##--> mlm.full is the best-fitting model and has all the weight

##
## some other model checking/diagnostics
##
# visualise r-hat scores
rhats = rhat(mlm.full)
color_scheme_set("brightblue")
mcmc_rhat(rhats)  # these all look pretty good

# visualise effect sample size ratios
ratios = neff_ratio(mlm.full)
mcmc_neff(ratios)

# check that none of the n_eff / N ratios are < 0.5
ratios[ ratios < 0.5 ]

## all looks good


#######################################################################
## Plot model coefficients - site-specific intercepts and slopes
##
draws = as.array(mlm.full)
draws = draws[,, grep( "Pop:Ego", dimnames(draws)$parameters, invert = T )]  # drop intercepts for egos nested in sites

# for each population, calculate its intercept and slopes
pops = unique(as.character(d.long$Pop))

# add intercepts for each site to population-average intercept
draws.int        = draws[, grep(paste( paste0("Intercept.*",    pops), collapse="|"), dimnames(draws)$parameters)] + draws[, "(Intercept)"]
# add slopes for each site to average slopes
draws.slope.r    = draws[, grep(paste( paste0("b\\[r Pop\\:",   pops), collapse="|"), dimnames(draws)$parameters)] + draws[, "r"]
draws.slope.hg   = draws[, grep(paste( paste0("b\\[Same.*",     pops), collapse="|"), dimnames(draws)$parameters)] + draws[, "SameHerdingGroup"]
draws.slope.r_hg = draws[, grep(paste( paste0("b\\[r\\:Same.*", pops), collapse="|"), dimnames(draws)$parameters)] + draws[, "r:SameHerdingGroup"]

# prettify names and calculate quantiles for plotting
dimnames(draws.int)$parameters        = paste0("Intercepts_",pops)
dimnames(draws.slope.r)$parameters    = paste0("Relatedness_",pops)
dimnames(draws.slope.hg)$parameters   = paste0("Herding group_",pops)
dimnames(draws.slope.r_hg)$parameters = paste0("Interaction_",pops)

draws.site = cbind(draws.int, draws.slope.r, draws.slope.hg, draws.slope.r_hg)

quants = apply(draws.site, 2, quantile, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))  # quantiles of posterior

# convert to data frame, rename the columns, turn the rownames into a 'parameter' column, add a grouping variable for intercepts/slopes
quants = as.data.frame(t(quants))
names(quants) = c("lwr.95", "lwr.80", "med", "upr.80", "upr.95")
quants = cbind(Parameter=rownames(quants), quants)  # convert row names into a 'Parameter' variable
rownames(quants) = NULL  # remove row names
quants = quants %>% 
  separate(Parameter, c("Type", "Site"), sep="_")

quants$Type = factor(quants$Type, levels=c("Intercepts", "Relatedness", "Herding group", "Interaction"))
quants$Site = factor( convert_pop(quants$Site) )

# min(quants$med)
# max(quants$med)

# plot estimates
# (use forcats::fct_rev() to get the site labels in alphabetical order from top to bottom)
ggplot(quants, aes(x=forcats::fct_rev(Site), y=med)) +
  geom_hline(yintercept=0, colour = "grey80", linetype=2) +  # line of no effect
  geom_linerange(aes(ymin=lwr.80, ymax=upr.80, colour=Site), size=1.5) +
  geom_linerange(aes(ymin=lwr.95, ymax=upr.95, colour=Site)) +
  geom_point(aes(colour=Site), fill="white", size=2, stroke=1.5, shape = 21) +
  
  facet_wrap(~Type, strip.position = "left", ncol=1) +  # add grouping label (source: http://stackoverflow.com/a/36337286)
  coord_flip() +
  
  scale_color_brewer(type="qual", palette="Dark2") +
  
  xlab("") +
  ylab("Log-odds") +
  
  theme_bw() +
  theme(panel.spacing = unit(0, "lines") 
        ,strip.background = element_blank()
        ,strip.placement = "outside"
        ,plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.background = element_blank()
        ) +
  theme(legend.position="none")

ggsave(filename=file.path(plots.dir, "gifts - site-specific coefficients.png"), width=150, height=200, units="mm")
ggsave(filename=file.path(plots.dir, "gifts - site-specific coefficients.pdf"), width=150, height=200, units="mm")


##
## make a table of all coefficients
##
quants_table = quants %>% 
  select(Site, Parameter=Type, med, lwr.95, upr.95) %>% 
  mutate(Parameter = str_replace(Parameter, "slopes", ""), Site = as.character(Site)) %>% 
  arrange(Site)

write_csv(quants_table, file.path(results.dir, "gifts - param estimates.csv"))


#######################################################################
## Variance and correlation components
## compare in intercept-only (clustering only) model to best model
##
VPCs = function(m) {
  # get variances for varying intercepts and slopes, and add on the population-average intercept variance
  vc = VarCorr(m) %>% 
    dplyr::as_data_frame() %>% 
    filter(is.na(var2)) %>%      # get rid of covariances
    select(-var2) %>%            #... and unneeded columns
    add_row(grp="Avg", var1="(Intercept)",  # add population-average intercept variance to the beginning of this dataframe
            vcov=m$covmat["(Intercept)", "(Intercept)"], .before=1)
  
  # calculation variance partition coefficients and report back
  vc %>% 
    mutate(vpc = vcov / sum(vcov))
}

# calculate VPCs for all models
vpc_list = lapply(mlms, function(m) VPCs(m))
names(vpc_list) = models_list

##
## model summaries
##
# delta LOO information criteria
IC.list <- abs( unlist( sapply(loos, function(x) x[ grep("^elpd", names(x)) ]) ) )
dIC <- IC.list - min( IC.list )
names(dIC) = sapply(loos, function(x) attr(x, "model"))

# numbers of groups and observations in each model
model_summaries = sapply(mlms, function(x) c(ngrps(x), N=nobs(x)))
model_summaries = rbind(model_summaries, dIC, mlm_weights)  # tag on delta ICs

##
## save everything to a single Excel workbook - one sheet for VPC table and one sheet containing model summaries
##
wb = createWorkbook()

sheet_vpc = createSheet(wb, "VPCs")
start_col = 1
for (v in vpc_list) {
  addDataFrame(as.data.frame(v), sheet=sheet_vpc, startColumn=start_col, row.names=F)
  start_col = start_col + ncol(v)
}

sheet_sum = createSheet(wb, "Summaries")
addDataFrame(model_summaries, sheet=sheet_sum, startColumn=1, row.names=F)

saveWorkbook(wb, file.path(results.dir, "Model summaries.xlsx"))

rm(dIC, IC.list, model_summaries)
rm(wb, sheet_vpc, start_col, v, sheet_sum)
rm(mlms)


#######################################################################
## Plot predicted probabilities of gifts
##
new.data = expand.grid(
  r = seq(from=0, to=0.5, length.out=30),
  SameHerdingGroup = 0:1,
  Pop = unique(as.character(d.long$Pop)),
  Ego.ID = 0, GiftGiven = 0
)

# ~ (... | Pop) ==> average over egos in sites
post = posterior_predict(mlm.full, new.data, re.form = ~ (r * SameHerdingGroup | Pop))
pred = posterior_linpred(mlm.full, new.data, re.form = ~ (r * SameHerdingGroup | Pop), transform=T)

quants = apply(post, 2, quantile, probs = c(0.025, 0.5, 0.975))  # quantiles over mcmc samples
quants2 = apply(pred, 2, quantile, probs = c(0.025, 0.5, 0.975))  # quantiles over mcmc samples
row.names(quants) = c("sim.lwr", "sim.med", "sim.upr")
row.names(quants2) = c("lwr", "GiftGiven.pred", "upr")

new.data = cbind(new.data, t(quants), t(quants2))

# rename and re-order study sites
new.data$PopName = factor( convert_pop(new.data$Pop) )

# dummy data to give Finnmark the same y limits and Tibet the same y limits
# plots points at r = 1 (which will be hidden by the x limits in scale_x_log10)
ylims = data.frame(GiftGiven.pred=c(0.9, 0.9, 0.3, 0.3, 0.3, 0.3), PopName=levels(new.data$PopName), 
                   r=1, SameHerdingGroup=0)

ggplot(new.data, aes(x=r, y=GiftGiven.pred)) +
  geom_line(aes(colour=factor(SameHerdingGroup))) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=factor(SameHerdingGroup)), alpha=0.2) +
  
  geom_blank(data=ylims) +  # plot dummy data to give facets the same y limits by country

  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  
  xlab("Coefficient of relatedness") +
  ylab("Probability of gift giving") +
  
  # plot relatedness on log scale with breaks at familiar points
  scale_x_log10(breaks=c(0.0001, 0.0039, 0.0078, 0.015, 0.031, 0.063, 0.125, 0.25, 0.5), limits=c(0.0115, 0.5)) +  

  facet_wrap(~PopName, scales="free_y", ncol=2) +
  common_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # rotate x-axis labels
  )
  
ggsave(filename=file.path(plots.dir, "Gifts - predicted - multilevel models.png"), width=175, height=200, units="mm")
ggsave(filename=file.path(plots.dir, "Gifts - predicted - multilevel models.pdf"), width=175, height=200, units="mm")
