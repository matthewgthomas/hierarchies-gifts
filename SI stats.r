source("init if necessary.r")

library(dplyr)

# number/proportion of Saami herders we have herd size data for
people %>% 
  filter(Pop %in% c("karasjok", "kautokeino")) %>% 
  mutate(HerdData = ifelse(is.na(HerdSize), 0, 1)) %>% 
  summarise(Data_for = sum(HerdData), N = n(), percent = (sum(HerdData) / n()) * 100)

# minimum herd sizes in each site
people %>% 
  group_by(Pop) %>% 
  filter(HerdSize > 0) %>% 
  summarise(min_herd = min(HerdSize))

# sex ratios by country
# we don't know the gender of one person in Kautokeino - they are coded as 'o' for 'other'
people %>% 
  mutate(Country = ifelse(str_sub(Pop, start=1, end=1) == "k", "Norway", "China")) %>% 
  group_by(Country, Sex) %>% 
  ftable(Sex ~ Country, data=.) %>% 
  prop.table(margin=1) %>% 
  round(4)

# sex ratios by site
people %>% 
  group_by(Pop, Sex) %>% 
  ftable(Sex ~ Pop, data=.) %>% 
  prop.table(margin=1) %>% 
  round(4)

# age histograms
people %>% 
  mutate(Country = ifelse(str_sub(Pop, start=1, end=1) == "k", "Norway", "China")) %>%
  ggplot(aes(x=Age)) +
  geom_histogram(binwidth=5) +
  facet_wrap(~Country, scales="free_y")
