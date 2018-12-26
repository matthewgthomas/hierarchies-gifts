##
## Plot descriptive stats
##
library(futile.logger)
library(xlsx)
library(tidyverse)
source("init if necessary.r")

flog.info("Plotting descriptive stats")


############################################################################################
## Gifts split by whether or not they were given to members of the same herding group:
##
dyads.subset %>% 
  filter(GiftGiven==1) %>% 
  
  ggplot(aes(x=factor(SameHerdingGroup))) +
    geom_histogram(stat="count") +
    
    scale_x_discrete(labels=c("0" = "No", "1" = "Yes")) +
    xlab("Member of same herding group?") +
    ylab("Number of gifts") +
    
    facet_wrap(~PopName) +
    common_theme
    
ggsave(filename=file.path(plots.dir, "gifts by herding group membership.png"), width=150, height=150, units="mm")
ggsave(filename=file.path(plots.dir, "gifts by herding group membership.pdf"), width=150, height=150, units="mm")


############################################################################################
## Ranges of gifts by site
##
people %>% 
  group_by(PopName) %>% 
  summarise(
    # summarise number of gifts received
    n.mean = mean(gift.deg.in), n.median = median(gift.deg.in),
    n.min  = min(gift.deg.in),  n.max = max(gift.deg.in),
    # summarise amounts received
    amount.mean = mean(gifts.total), amount.median = median(gifts.total),
    amount.min  = min(gifts.total),  amount.max = max(gifts.total)) %>% 
  write_csv(path=file.path(results.dir, "gifts summary.csv"))


############################################################################################
## Kin and non-kin in the same herding group and other groups:
##
dyads.r = dyads.subset %>% 
  mutate(r.bin = cut(r, c(0, 0.0039, 0.0078, 0.015, 0.031, 0.063, 0.125, 0.25, 0.5, 1), right=F))

dyads.r$SameHerdingGroupText = factor(dyads.r$SameHerdingGroup)
levels(dyads.r$SameHerdingGroupText) = c("No", "Yes")

ggplot(dyads.r, aes(x=r.bin)) +
  geom_histogram(aes(fill=SameHerdingGroupText), position="dodge", stat="count", alpha=0.75) +
  
  scale_y_log10(breaks=c(1, 10, 100, 1000, 10000)) +
  
  scale_fill_manual(values = cols) +
  xlab("Coefficient of relatedness") +
  ylab("Frequency") +
  
  guides(fill= guide_legend("Same herding group?")) +
  facet_wrap(~PopName, scales="free_y") +
  
  common_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # rotate x-axis labels
   ,legend.position="bottom"
  )

ggsave(filename=file.path(plots.dir, "relatedness to herding group.png"), width=250, height=150, units="mm")
ggsave(filename=file.path(plots.dir, "relatedness to herding group.pdf"), width=250, height=150, units="mm")

rm(dyads.r)


############################################################################################
## Proportion of gifts going to a close relative (r>=0.5), any other relative, non-kin, split by population
##
# divide into relatedness categories
dyads.gifts = dyads.subset %>% 
  filter(GiftGiven == 1) %>% 
  mutate(r.desc = factor( case_when(
    .$r >= 0.5              ~ "Immediate family",
    .$r >= 0.25 & .$r < 0.5 ~ "Other close kin",
    .$r > 0 & .$r < 0.25    ~ "Distant kin",
    .$r == 0                ~ "Non-kin"
  ), levels=c("Non-kin", "Distant kin", "Other close kin", "Immediate family"))) %>% 
  
  # calculate proportions of gifts given to different classes of kin
  group_by(PopName, r.desc) %>% 
  summarise(n = n()) %>% 
  mutate(Prop = n / sum(n))

ggplot(dyads.gifts, aes(x=r.desc, y=Prop)) +
  geom_bar(stat="identity") +
  
  xlab("") +
  ylab("Proportion of gifts") +
  facet_wrap(~PopName) +
  
  common_theme +
  theme( axis.text.x = element_text(angle = 45, hjust = 1) )  # rotate x-axis labels

ggsave(filename=file.path(plots.dir, "gifts by relatedness.png"), width=150, height=150, units="mm")
ggsave(filename=file.path(plots.dir, "gifts by relatedness.pdf"), width=150, height=150, units="mm")

rm(dyads.gifts)


############################################################################################
## Distribution of herd sizes for everyone in the sample
##
# get median herd sizes in each site and draw on graph
median.cattle = people %>% 
  group_by(PopName) %>% 
  summarise(Median.cattle = median(HerdSize, na.rm=T))

ggplot(people, aes(x=HerdSize)) +
  geom_histogram(binwidth = 50, alpha=0.5) +
  geom_vline(aes(xintercept = Median.cattle), median.cattle) +
  
  xlab("Number of livestock") +
  ylab("Frequency") +
  
  facet_wrap(~PopName, scales = "free") +
  common_theme

ggsave(filename=file.path(plots.dir, "distribution of herd sizes.png"), width=200, height=150, units="mm")
ggsave(filename=file.path(plots.dir, "distribution of herd sizes.pdf"), width=200, height=150, units="mm")

rm(median.cattle)


############################################################################################
## summary stats table
##
# calculate each person's relatedness to their herding group (in each site)
group.r = dyads %>% 
  filter(SameHerdingGroup == 1) %>%        # only want within-group dyads
  unite(Pop.HG, PopName, Ego.HG) %>%       # group IDs might be repeated across sites, so make a unique identifier
  select(DyadID, Pop.HG, r) %>%            # keep only necessary variables...
  distinct() %>%                           #... and remove repeated dyads (essentially making an undirected network)
  group_by(Pop.HG) %>%                     # summary stats for each site/herding group
  summarise(r.mean   = mean(r, na.rm=T),
            r.median = median(r, na.rm=T),
            r.sd     = sd(r, na.rm=T),
            r.se     = sd(r) / sqrt(length(r)),
            N.dyads  = length(r))
# separate(Pop.HG, c("PopName", "HG"), sep="_")  # split the Pop.HG unique identified (step 2) back into its separate components

# how many people in each herding group (that we interviewed)?
group.n = people %>% 
  unite(Pop.HG, PopName, HG) %>%
  group_by(Pop.HG) %>%
  summarise(N = n())

# merge group sizes into group relatedness table
group.r = group.r %>% 
  full_join(group.n, by="Pop.HG") %>%                 # use full_join because there are some groups in group.n who we don't have relatedness info for (because there are no within-group dyads)
  separate(Pop.HG, c("PopName", "HG"), sep="_") %>%   # split the Pop.HG unique identified back into its separate components
  arrange(PopName, N)  # order by number of members

# grand mean of relatedness within each site (and some other summary stats)
group.r.sum = group.r %>% 
  group_by(PopName) %>% 
  summarise(n.groups = n(),
            mean.N = mean(N, na.rm=T),
            sd.N = sd(N, na.rm=T),
            grand.mean.r = mean(r.mean, na.rm=T))

# save this summary
write.csv(group.r.sum, file=file.path(results.dir, "herding groups summary.csv"), row.names=F)

# summarise total N
sample.summary = dyads %>% 
  group_by(PopName) %>% 
  summarise(n = length(unique(Ego.ID)))

# summarise people who played gift games
gifts.summary = dyads %>% 
  group_by(PopName) %>% 
  filter(GiftGiven==1) %>% 
  summarise(n.givers = length(unique(Ego.ID)), n.gifts = sum(GiftGiven))

# summarise ages and genders
age.summary = people %>% 
  select(HerderID, PopName, Age) %>% 
  na.omit %>% 
  group_by(PopName) %>% 
  summarise(age.mean = mean(Age),
            age.sd   = sd(Age))

sex.summary = people %>% 
  select(HerderID, PopName, Sex) %>% 
  na.omit %>% 
  count(PopName, Sex) %>% 
  spread(Sex, n, fill=0)
  
# summarise herd sizes
herds.summary = people %>% 
  select(HerderID, PopName, HerdSize, cattle.z, Age.z, Sexb, gift.deg.in, gift.betweenness, gift.eigen) %>% 
  na.omit %>% 
  group_by(PopName) %>% 
  summarise(herd.mean = mean(HerdSize),
            herd.sd = sd(HerdSize))

# merge in gifts and herding groups summaries
sample.summary = sample.summary %>% 
  left_join(gifts.summary, by="PopName") %>% 
  left_join(age.summary,   by="PopName") %>% 
  left_join(sex.summary,   by="PopName") %>% 
  left_join(group.r.sum,   by="PopName") %>% 
  left_join(herds.summary, by="PopName")

# calculate totals
sample.totals = sample.summary %>% 
  summarise(PopName="Totals", n=sum(n), n.givers=sum(n.givers), n.gifts=sum(n.gifts), 
            age.mean=99, age.sd=99, f=sum(f), m=sum(m), o=sum(o),
            n.groups=sum(n.groups), mean.N=-99, sd.N=-99, grand.mean.r=-99, herd.mean=-99, herd.sd=-99)  # treat this as missing -- don't want to total the averages -- but NA doesn't work properly when opening in Excel

# save summary
wb = createWorkbook()
sheet_sample = createSheet(wb, "sample descriptives")
addDataFrame( as.data.frame(bind_rows(sample.summary, sample.totals)), sheet=sheet_sample, row.names=F)
saveWorkbook(wb, file.path(results.dir, "Sample summaries.xlsx"))
