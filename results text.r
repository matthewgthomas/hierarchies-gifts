#' ---
#' title: Results
#' output: word_document
#' ---

#+ echo=F
library(dplyr)
library(xlsx)

source("init if necessary.r")

n.herders = length(unique(dyads$Ego.ID))
n.givers  = length(unique(dyads.subset$Ego.ID))

finnmark_gifts.total = dyads %>% 
  filter(grepl("Finnmark*", PopName)) %>% 
  group_by(SameHerdingGroup) %>% 
  summarise(gifts = sum(GiftGiven)) %>% 
  mutate(pc = round( (gifts / sum(gifts)) * 100, 1))

finnmark_gifts.to_group = finnmark_gifts.total %>% filter(SameHerdingGroup==1) %>% select(pc) %>% as.numeric(.)

tibet_gifts.total = dyads %>% 
  filter(grepl("Tibet*", PopName)) %>% 
  group_by(SameHerdingGroup) %>% 
  summarise(gifts = sum(GiftGiven)) %>% 
  mutate(pc = round( (gifts / sum(gifts)) * 100, 1))

tibet_gifts.to_group = tibet_gifts.total %>% filter(SameHerdingGroup==1) %>% select(pc) %>% as.numeric(.)

# reciprocities
finnmark.stats = net.stats %>% 
  filter(grepl("^k.*", Site)) %>% 
  summarise(recip.min = round( min(reciprocity)*100, 2), recip.max = round( max(reciprocity)*100, 2),
            assort.min = round(min(assortativity.hg), 2), assort.max = round(max(assortativity.hg), 2))

tibet.stats = net.stats %>% 
  filter(!grepl("^k.*", Site)) %>% 
  summarise(recip.min = round( min(reciprocity)*100, 2), recip.max = round( max(reciprocity)*100, 2),
            assort.min = round(min(assortativity.hg), 2), assort.max = round(max(assortativity.hg), 2))

#' In total, `r n.givers` herders living in six study sites gave `r sum(d.long$GiftGiven)` gifts. Overall our sample contained `r n.herders` herders in `r nrow(dyads.subset)` dyads. Reindeer herders in Finnmark gave `r finnmark_gifts.to_group`% of their `r sum(finnmark_gifts.total$gifts)` gifts to members of the same herding group, while yak herders in Tibet gave `r tibet_gifts.to_group`% of `r sum(tibet_gifts.total$gifts)` gifts to members of the same group (Figure...).
#'
#' Finnmark: `r finnmark.stats$recip.min`-`r finnmark.stats$recip.max`% of gifts were reciprocated. Assortment on herding group: `r finnmark.stats$assort.min`-`r finnmark.stats$assort.max`

#+ echo=F
post = as.data.frame(m.herd.deg)
post.sex = quantile(post$Sexb, probs=c(0.5, 0.025, 0.975))
post.age = quantile(post$Age.z, probs=c(0.5, 0.025, 0.975))

#' Females had smaller herds compared to males ($\Beta$=`r round(post.sex[1], 3)`; 95% credible interval [`r round(post.sex[2], 3)`, `r round(post.sex[3], 3)`]), while older herders owned more animals ($\Beta$=`r round(post.age[1], 3)`; 95% CI [`r round(post.age[2], 3)`, `r round(post.age[3], 3)`]).

#+ echo=F
mod_kara = net.stats %>% filter(Site=="karasjok") %>% select(modularity) %>% as.numeric(.)
mod_kara_prop = net.stats.rnd %>% filter(Site=="karasjok" & modularity < mod_kara) %>% summarise(p = n() / nrow(net.stats.rnd)) %>% as.numeric(.)

mod_kauto = net.stats %>% filter(Site=="kautokeino") %>% select(modularity) %>% as.numeric(.)
mod_kauto_prop = net.stats.rnd %>% filter(Site=="kautokeino" & modularity < mod_kauto) %>% summarise(p = n() / nrow(net.stats.rnd)) %>% as.numeric(.)

#' `r round(mod_kara_prop * 100, 1)`% of random modularity scores < observed score in Karasjok; `r round(mod_kauto_prop * 100, 1)`% of random modularity scores < observed score in Kautokeino

#+ echo=F
vpc.herd.sum = round(vpc.herd.deg$vpc * 100, 1)

#' Between-subject differences accounted for almost all variance in predicting standardised herd size (`r vpc.herd.sum[2]`%); there was almost no variation between sites (`r vpc.herd.sum[1]`%).

#+ echo=F
# calculate proportions of gifts given to non-kin and kin
dyads.gifts = dyads.subset %>% 
  filter(GiftGiven == 1) %>% 
  mutate(r.desc = factor( case_when(
    .$r > 0  ~ "Kin",
    .$r == 0 ~ "Non-kin"
  ), levels=c("Non-kin", "Kin"))) %>% 
  group_by(PopName, r.desc) %>% 
  summarise(n = n()) %>% 
  mutate(Prop = n / sum(n)) %>% 
  ungroup()

kauto_prop = dyads.gifts %>%
  filter(PopName=="Norway: Kautokeino" & r.desc=="Kin") %>% 
  select(Prop) %>% 
  as.numeric(.)

kara_prop = dyads.gifts %>%
  filter(PopName=="Norway: Karasjok" & r.desc=="Non-kin") %>% 
  select(Prop) %>% 
  as.numeric(.)

tibet_props = dyads.gifts %>%
  filter(grepl("Tibet", PopName) & r.desc=="Non-kin") %>% 
  select(Prop) %>% 
  range(.)

#' Proportion of gifts to kin in Kautokeino: `r round(kauto_prop * 100, 1)`%
#' 
#' Proportion of gifts to non-kin in Karasjok: `r round(kara_prop * 100, 1)`%
#' 
#' Range of non-kin gift proportions in Tibet: `r round(tibet_props[1] * 100, 1)`% &ndash; `r round(tibet_props[2] * 100, 1)`%

#+ echo=F
library(ineq)

ginis = dyads %>% 
  select(Alter.ID, Alter.HerdSize, PopName) %>% 
  distinct() %>% 
  group_by(PopName) %>% 
  summarise(cv = sd(Alter.HerdSize, na.rm=T) / mean(Alter.HerdSize, na.rm=T),
            gini = ineq(Alter.HerdSize, type="Gini"))

gini_finnmark = ginis %>% 
  filter(grepl("Finnmark", PopName)) %>% 
  select(gini) %>% 
  range(.)

gini_tibet = ginis %>% 
  filter(grepl("Tibet", PopName)) %>% 
  select(gini) %>% 
  range(.)

#' Finnmark Gini range: `r round(gini_finnmark[1], 3)` &ndash; `r round(gini_finnmark[2], 3)`
#' 
#' Tibet Gini range: `r round(gini_tibet[1], 3)` &ndash; `r round(gini_tibet[2], 3)`

