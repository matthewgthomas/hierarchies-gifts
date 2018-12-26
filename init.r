##
## load dyadic and individual data
##
library(futile.logger)
library(data.table)
library(tidyverse)
library(stringr)
library(igraph)

# decide whether to save .csv files containing individuals and dyads
# (if you wanted to load data directly from these files rather than recalculating every time by running this code)
save_data_files = F

# some folders
data.dir = "data"
results.dir = "results"
plots.dir = "plots"
models.dir = "models"
docs.dir = "docs"  # for the map

# create the directories if they don't exist
if (!dir.exists(plots.dir))    # plots
  dir.create(plots.dir)
if (!dir.exists(results.dir))  # results
  dir.create(results.dir)
if (!dir.exists(models.dir))   # models
  dir.create(models.dir)
if (!dir.exists(docs.dir))     # map
  dir.create(docs.dir)

# some filenames that are used in more than one source file
gift_models_file = file.path(models.dir, "gift analyses - multilevel models.RData")
comparison_data_file = file.path(models.dir, "gift analyses - model comparison.RData")

flog.info("Preparing data...")


############################################################################
## helper function to convert site names in the data to something more readable
##
convert_pop = function(p)
{
  dplyr::case_when(
    p == "cairima"    ~ "Tibet: Cairima",
    p == "doulong"    ~ "Tibet: Doulong",
    p == "jilehe"     ~ "Tibet: Jilehe",
    p == "tawa"       ~ "Tibet: Tawa",
    p == "karasjok"   ~ "Finnmark: Karasjok",
    p == "kautokeino" ~ "Finnmark: Kautokeino"
  )
}


############################################################################
## load nodes and edges
##
node.files  = list.files(data.dir, pattern="nodes.*csv", full.names=F)
r.files     = list.files(data.dir, pattern="edges - r.*csv", full.names=F)
gifts.files = list.files(data.dir, pattern="edges - gifts.*csv", full.names=F)

covars = c("HG", "HerdSize", "Age", "Sex")

# x = 1
# loop over all the files (in alphabetical order)
for (x in 1:length(node.files))
{
  herders   = as.data.table( read.csv(file.path(data.dir, node.files[x]), stringsAsFactors = F) )
  herders.r = as.data.table( read.csv(file.path(data.dir, r.files[x]), stringsAsFactors = F) )
  gifts     = as.data.table( read.csv(file.path(data.dir, gifts.files[x]), stringsAsFactors = F) )
  
  site = substring(node.files[x], 9, nchar(node.files[x]) - 4)  # "nodes - " is 8 characters so start at 9th, ".csv" is final 4
  
  setkey(herders,   HerderID)
  setkey(herders.r, Ego, Alter)
  setkey(gifts,     Ego, Alter)
  
  # recode sex if it's not already a character -- (1 is female)
  if (is.numeric(herders$Sex))
    herders[, Sex := ifelse(Sex == 1, "f", "m")]

  
  #########################################################################
  ## Make main dyadic data frame in wide format
  ##
  dyads = data.table(expand.grid(Ego=herders$HerderID, Alter=herders$HerderID))
  setkeyv(dyads, c("Ego", "Alter"))
  # remove dyads where ego==alter
  dyads = dyads[Ego != Alter]
  
  dyads$Pop = site
  
  ##
  ## Set dyad IDs for each ego/alter pair
  ##
  # set dyad ID to be the smallest of ego/alter followed by the largest
  # dyads[, DyadID := ifelse(Ego < Alter, paste(Ego, Alter, sep=""), paste(Alter, Ego, sep=""))]
  # dyads[, DyadID := as.integer(DyadID)]  # paste() makes it character; convert to number
  
  ##
  ## Merge in gifts
  ##
  dyads = gifts[dyads]
  setkeyv(dyads, c("Ego", "Alter"))  # need to reset keys after merge
  
  # remove NAs
  dyads[is.na(Amount), Amount := 0]
  
  # create binary response variable for gifts
  dyads[, GiftGiven := 0]
  dyads[Amount > 0, GiftGiven := 1]
  
  ##
  ## Merge in relatedness
  ##
  # merge coefficients of relatedness
  dyads = herders.r[dyads]
  # set NAs to zero
  dyads[is.na(r), r := 0]
  
  ##
  ## Merge in covariates for ego and alter
  ##
  setkey(dyads, Ego)
  dyads = herders[dyads]
  setnames(dyads, "HerderID", "Ego")  # sort out column names
  
  # repeat for alter
  setkey(dyads, Alter)
  dyads = herders[dyads]
  setnames(dyads, "HerderID", "Alter")  # sort out column names
  
  # sort out covariates' names
  setnames(dyads, covars, paste("Alter", covars, sep="."))
  setnames(dyads, paste("i", covars, sep="."), paste("Ego", covars, sep="."))
  
  # setkey(dyads, Ego, Alter)
  
  ##
  ## same herding group?
  ##
  dyads[, SameHerdingGroup := 0]
  dyads[Ego.HG==Alter.HG, SameHerdingGroup := 1]  # belong to same herding group?
  
  ##
  ## ensure against duplicate ID numbers
  ##
  # setkey(dyads, DyadID)
  # dyads[, DyadID := paste0(substr(Pop, 1, 3), DyadID)]
  
  # some ego/alter ID numbers are duplicated across sites; prepend the first two letters of the site name to make IDs properly unique
  dyads = dyads %>% 
    mutate(tmpPop = substr(Pop, 1, 3)) %>%                  # get first two letters of site
    unite(Ego.ID, tmpPop, Ego, sep="", remove=F) %>%        # prepend site to ego ID
    select(-Ego) %>%                                        # get rid of original 'ego' column
    unite(Alter.ID, tmpPop, Alter, sep="", remove=T) %>%    # prepend site to alter ID
    as.data.table()
  
  
  ##
  ## Set dyad IDs for each ego/alter pair
  ## - using 'str_extract()' to help put the alphanumeric IDs in numerical order
  ##
  dyads[, DyadID := ifelse( as.integer(str_extract(Ego.ID, "[0-9]+")) < as.integer(str_extract(Alter.ID, "[0-9]+")),
                            paste0(Ego.ID, Alter.ID), paste0(Alter.ID, Ego.ID))]
  
  setkey(dyads, DyadID)
  
  
  ##
  ## Sort out column order and save dyadic data file
  ##
  # reorder the columns into something more sensible
  setcolorder(dyads, c("DyadID", 
                       "Ego.ID", paste("Ego", covars, sep="."),
                       "Alter.ID", paste("Alter", covars, sep="."),
                       "r", "SameHerdingGroup", "GiftGiven", "Amount",
                       "Pop"
  ))
  # setkey(dyads, DyadID)
  # names(dyads)
  
  assign( paste0("dyads.", site), dyads )
  if (save_data_files) write_csv(dyads, file.path(data.dir, paste0("dyads - ", site, ".csv")))

  
  #########################################################################
  ## herders
  ##
  herders$Pop = site
  
  ##
  ## Calculate centrality in gift network
  ##
  # create gift network
  gifts = subset(gifts, Ego   %in% herders$HerderID & Alter %in% herders$HerderID)
  
  g.gifts = graph.data.frame(gifts, vertices=herders, directed=T)
  
  # in-degree (number of gifts received)
  deg.in  = degree(g.gifts, mode="in")
  herders$gift.deg.in = deg.in[ as.character(herders$HerderID) ]
  
  # total amount of currency received
  amounts = gifts %>% 
    group_by(HerderID = Alter) %>% 
    summarise(gifts.total = sum(Amount))
  
  herders = herders %>% 
    left_join(amounts, by="HerderID")
  
  herders = mutate(herders, gifts.total = ifelse(is.na(gifts.total), 0, gifts.total))  # remove NAs
  
  # eigenvector centrality and betweenness
  herders$gift.eigen = eigen_centrality(g.gifts, directed=T, scale=F, weights=NA)$vector
  herders$gift.betweenness = betweenness(g.gifts, directed=T, normalized=F, weights=NA)
  
  ##
  ## ensure against duplicated ID numbers
  ##
  # some ego/alter ID numbers are duplicated across sites; prepend the first two letters of the site name to make IDs properly unique
  herders = herders %>% 
    mutate(tmpPop = substr(Pop, 1, 3)) %>%     # get first two letters of site
    unite(HerderID, tmpPop, HerderID, sep="", remove=T)  # prepend site to ego ID
 
  assign( paste0("herders.", site), herders )
  
  flog.info(paste0("Finished ", site))
}

# clean up
rm(node.files, gifts.files, r.files)
rm(dyads, gifts, herders, herders.r)
rm(covars, deg.in, amounts, g.gifts, site, x)


#########################################################################
## merge all the separate dyads and herders data frames into one
##
d_vars = grep("dyads\\.", ls(), value=T)

dyads = rbindlist(mget( d_vars ))

# variable formatting
dyads[, SameHerdingGroup := as.integer(SameHerdingGroup)]
dyads[, GiftGiven := as.integer(GiftGiven)]
dyads[, Pop := as.factor(Pop)]

# there are a couple of people whose sex we don't know; mark them as 'o' rather than NA
dyads[, Ego.Sex   := ifelse(is.na(Ego.Sex),   "o", Ego.Sex)]
dyads[, Alter.Sex := ifelse(is.na(Alter.Sex), "o", Alter.Sex)]
dyads[, Ego.Sex   := as.factor(Ego.Sex)]
dyads[, Alter.Sex := as.factor(Alter.Sex)]

dyads[, PopName := convert_pop(Pop)]  # add friendly name for each site

rm(list = d_vars)
rm(d_vars)


############################################################################
## Subsets of dyadic data
##
# keep only Egos who played the gift game in each site
gift_givers = dyads %>% 
  filter(GiftGiven==1) %>% 
  select(Ego.ID) %>% 
  distinct()

dyads.subset = dyads %>% 
  filter(Ego.ID %in% gift_givers$Ego.ID)

# dataframe for multilevel models
d.long = subset(dyads.subset, select=c(Ego.ID, Alter.ID, Pop, r, SameHerdingGroup, GiftGiven))
d.long = na.omit(d.long)


############################################################################
## Create people table
##
h_vars = grep("herders\\.", ls(), value=T)
people = rbindlist(mget( h_vars ))
rm(list = h_vars)
rm(h_vars)

# sort out variable formats
people[, Pop := as.factor(Pop)]
people[, Sex := as.factor(Sex)]

people[, PopName := convert_pop(Pop)]  # add friendly name for each site

# standardise number of cattle (grouped within sites)
people[, cattle.z := scale(HerdSize), by="Pop"]
people[, Age.z := scale(Age), by="Pop"]
people[, Sexb := ifelse(Sex=="m", 0, 1)]


############################################################################
## Save
##
if (save_data_files) {
  write_csv(dyads,  file.path(data.dir, "dyads.csv"))
  write_csv(people, file.path(data.dir, "people.csv"))
  
  gifts = dyads %>% 
    filter(GiftGiven > 0) %>% 
    select(Ego = Ego.ID, Alter = Alter.ID, weight = Amount)
  
  g = graph.data.frame(gifts, vertices = people)
  # plot(g)
  write_graph(g, file = file.path(data.dir, "gifts.graphml"), format = "graphml")
}
