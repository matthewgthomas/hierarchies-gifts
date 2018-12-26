##
## social network analyses of gifts
##
library(futile.logger)
library(igraph)
library(tidyverse)

source("init if necessary.r")


####################################################################
## Prepare data
##
flog.info("SNA: preparing data")

# make a list of herder nodes
nodes = dyads %>% 
  select(starts_with("Ego"), Pop, PopName) %>% 
  distinct()

names(nodes) = gsub("Ego\\.", "", names(nodes))  # get rid of "Ego." at the start of column names

# sort out column types
#nodes$Sex = factor(nodes$Sex)
#nodes$HG = factor(nodes$HG)

# make an edge list of gifts
gifts = dyads %>% 
  filter(GiftGiven==1) %>% 
  select(Ego.ID, Alter.ID, Amount)
gifts$weight = gifts$Amount

g.gifts = graph.data.frame(gifts, vertices=nodes, directed=T)

study_sites = unique(nodes$Pop)

# no. gifts received
deg.in  = degree(g.gifts, mode="in")
nodes$NumGifts = deg.in[ as.character(nodes$ID) ]


####################################################################
## For each site, calculate various social network stats for observed and random nets
##
flog.info("SNA: simulating networks and calculating statistics...")

# set up data.frame for observed network stats
net.stats = data.frame(Site=rep("", length(study_sites)), n.obs=0,
                       clusters=0, reciprocity=0, density=0, transitivity=0, modularity=0,
                       assortativity.degree=0, assortativity.hg=0, assortativity.sex=0,
                       stringsAsFactors=F)

# set up data.frame for random nets
n.sims = 1000

net.stats.rnd = data.frame(Site=rep("", n.sims), sim=0, n.obs=0,
                           clusters=0, reciprocity=0, density=0, transitivity=0, modularity=0,
                           assortativity.degree=0, assortativity.hg=0, assortativity.sex=0,
                           stringsAsFactors=F)

i=1  # track which row of net.stats to write to
j=1  # track which row of net.stats.rnd to write to

for (pop in study_sites)
{
  #pop = unique(nodes$Pop)[1]
  
  # subset graph for this site only
  g.gifts.sub = induced_subgraph(graph=g.gifts, which(V(g.gifts)$Pop==pop))
  # get individuals in this sub-graph
  nodes.sub = subset(nodes, ID %in% V(g.gifts.sub)$name)
  
  # calc stats for this network
  net.stats[i,] = c(as.character(pop), 
                    length(E(g.gifts.sub)),
                    no.clusters(g.gifts.sub),
                    reciprocity(g.gifts.sub),
                    graph.density(g.gifts.sub),
                    transitivity(g.gifts.sub),
                    modularity(cluster_edge_betweenness(g.gifts.sub)),
                    assortativity(g.gifts.sub, types1 = nodes.sub$NumGifts),  # assorting on gifts?
                    0, #assortativity.nominal(g.gifts.sub, types=herding_groups, directed=T),  # problems with Tibet networks, so calculate separately below
                    #assortativity(g.gifts.sub, types1 = nodes.sub$cattle, types2 = nodes.sub$cattle),  # assorting on herd sizes? -- too many NAs
                    assortativity.nominal(g.gifts.sub, types=nodes.sub$Sex, directed=T)  # assorting on sex?
                    )
  
  ##
  ## we don't know herding group membership for some people in some of the Tibetan gift networks
  ## this causes problems (R crashes) when trying to calculated assortativity on herding group membership
  ## so, as a way to get around this, drop these people from the networks, just for calculating assortativity.nominal
  ##
  #which(!is.na(V(g.gifts.sub)$HG))
  g.gifts.sub.sub = induced_subgraph(graph=g.gifts.sub, which( !is.na(V(g.gifts.sub)$HG) ))
  nodes.sub.sub = subset(nodes, ID %in% V(g.gifts.sub.sub)$name)
  herding_groups = as.integer( factor( V(g.gifts.sub.sub)$HG ) )
  
  net.stats[i,]$assortativity.hg = assortativity.nominal(g.gifts.sub.sub, types=herding_groups, directed=T)
  
  ##
  ## calculate network stats for random networks with the same number of nodes and edges as in this site
  ##
  for (sim in 1:n.sims)
  {
    # make a random network
    g.rnd = sample_gnm(n=nrow(nodes.sub), m=length(E(g.gifts.sub)), directed=T)
    herding_groups = sample(1:length(unique(V(g.gifts.sub)$HG)), size = nrow(nodes.sub), replace=T)  # assign herding group membership randomly
    # sexes = sample(1:2, size = nrow(nodes.sub), replace=T)  # assign sex randomly
    # herding_groups[ floor(runif(20, 1, nrow(nodes.sub))) ] = NA
    
    net.stats.rnd[j,] = c(as.character(pop), sim,
                          length(E(g.rnd)),   # no. observations
                          no.clusters(g.rnd),
                          reciprocity(g.rnd),
                          graph.density(g.rnd),
                          transitivity(g.rnd),
                          modularity(cluster_edge_betweenness(g.rnd)),
                          assortativity.degree(g.rnd, directed=T),  # assorting on gifts?
                          assortativity.nominal(g.rnd, types=herding_groups, directed=T),  # assorting on herding groups?
                          assortativity.nominal(g.rnd, types=nodes.sub$Sex, directed=T))  # assorting on sex?
    j=j+1
  }
  
  i=i+1
  flog.info(paste0("Finished ", pop))
}

rm(g.gifts.sub, g.gifts.sub.sub, nodes.sub, nodes.sub.sub, i, herding_groups)
rm(g.rnd, j)

# clean up and save obseved stats
net.stats$n.obs = as.integer(net.stats$n.obs)
net.stats$clusters = as.integer(net.stats$clusters)
net.stats$reciprocity = as.numeric(net.stats$reciprocity)
net.stats$density = as.numeric(net.stats$density)
net.stats$transitivity = as.numeric(net.stats$transitivity)
net.stats$assortativity.degree = as.numeric(net.stats$assortativity.degree)
net.stats$assortativity.hg = as.numeric(net.stats$assortativity.hg)
net.stats$assortativity.sex = as.numeric(net.stats$assortativity.sex)
net.stats$modularity = as.numeric(net.stats$modularity)
# str(net.stats)

write.csv(net.stats, file=file.path(results.dir, "gift network statistics - observed.csv"), row.names=F)

# clean up random stats
net.stats.rnd$sim = as.integer(net.stats.rnd$sim)
net.stats.rnd$n.obs = as.integer(net.stats.rnd$n.obs)
net.stats.rnd$clusters = as.integer(net.stats.rnd$clusters)
net.stats.rnd$reciprocity = as.numeric(net.stats.rnd$reciprocity)
net.stats.rnd$density = as.numeric(net.stats.rnd$density)
net.stats.rnd$transitivity = as.numeric(net.stats.rnd$transitivity)
net.stats.rnd$assortativity.degree = as.numeric(net.stats.rnd$assortativity.degree)
net.stats.rnd$assortativity.hg = as.numeric(net.stats.rnd$assortativity.hg)
net.stats.rnd$assortativity.sex = as.numeric(net.stats.rnd$assortativity.sex)
net.stats.rnd$modularity = as.numeric(net.stats.rnd$modularity)
# str(net.stats.rnd)

write.csv(net.stats.rnd, file=file.path(results.dir, "gift network statistics - random networks.csv"), row.names=F)


####################################################################
## Plot observed and randomised network stats
##
flog.info("SNA: plotting network stats")

#' Function to plot network statistics (observed and simulated)
#'
#' @param net.stats.melt 
#' @param net.stats.rnd.melt 
#' @param cols, Number of columns of facets in the plot (default = 3)
#'
#' @return
#'
plot_net_stats = function(net.stats.melt, net.stats.rnd.melt, cols=3) {
  # format site/country names
  net.stats.melt$Site = stringr::str_to_title(net.stats.melt$Site)
  net.stats.melt$Country = ifelse( stringr::str_sub(net.stats.melt$Site, 1, 1) == "K", "Finnmark", "Tibet")
  
  net.stats.rnd.melt$Site = stringr::str_to_title(net.stats.rnd.melt$Site)
  net.stats.rnd.melt$Country = ifelse( stringr::str_sub(net.stats.rnd.melt$Site, 1, 1) == "K", "Finnmark", "Tibet")
  
  net.stats.melt$Site = paste0(net.stats.melt$Country, ": ", net.stats.melt$Site)
  net.stats.rnd.melt$Site = paste0(net.stats.rnd.melt$Country, ": ", net.stats.rnd.melt$Site)
  
  # plot them
  p = ggplot(net.stats.rnd.melt, aes(x=factor(Site), y=Value)) + 
    facet_wrap(~Measure, scales="free_y", ncol=cols) +
    geom_boxplot() + 
    geom_point(data=net.stats.melt, aes(x=factor(Site), y=Value), fill="red", shape=23, size=2) + 
    
    geom_blank(data=net.stats.melt) +
    
    ylab("Simulated/observed network statistics") +
    xlab("") +
    
    theme_bw() +
    #eliminates baground, gridlines, and chart border
    theme(
      plot.background = element_blank()
      ,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      ,panel.background = element_blank()
      ,axis.text=element_text(size=12)
      ,axis.title=element_text(size=12)
      ,axis.text.x = element_text(angle = 45, hjust = 1)  # rotate x-axis labels
    )
  p
}

##
## plot all measures
##
# prepare observed stats
net.stats.melt = net.stats %>% 
  select(Site, reciprocity, modularity, assortativity.hg, assortativity.sex) %>% 
  gather(Measure, Value, reciprocity:assortativity.sex, factor_key=T)

levels(net.stats.melt$Measure) = c("Reciprocity", "Modularity", "Herding group assortment", "Sex assortment")
net.stats.melt$Value = as.numeric(net.stats.melt$Value)

# prepare random stats
net.stats.rnd.melt = net.stats.rnd %>% 
  select(Site, reciprocity, modularity, assortativity.hg, assortativity.sex) %>% 
  gather(Measure, Value, reciprocity:assortativity.sex, factor_key=T)

levels(net.stats.rnd.melt$Measure) = c("Reciprocity", "Modularity", "Herding group assortment", "Sex assortment")

plot_net_stats(net.stats.melt, net.stats.rnd.melt, cols = 2)
ggsave(filename=file.path("plots", "network stats - observed and simulated - all.png"), width=250, height=200, units="mm")
ggsave(filename=file.path("plots", "network stats - observed and simulated - all.pdf"), width=250, height=200, units="mm")

##
## plot reciprociate and assorment on sex/herding groups only
##
# net.stats.melt = net.stats %>% 
#   select(Site, reciprocity, assortativity.hg:assortativity.sex) %>% 
#   gather(Measure, Value, reciprocity:assortativity.sex, factor_key=T)
# 
# levels(net.stats.melt$Measure) = c("Reciprocity", "Herding group assortment", "Sex assortment")
# #net.stats.melt$Value = as.numeric(net.stats.melt$Value)
# 
# # prepare random stats
# net.stats.rnd.melt = net.stats.rnd %>% 
#   select(Site, reciprocity, assortativity.hg:assortativity.sex) %>% 
#   gather(Measure, Value, reciprocity:assortativity.sex, factor_key=T)
# 
# levels(net.stats.rnd.melt$Measure) = c("Reciprocity", "Herding group assortment", "Sex assortment")
# 
# plot_net_stats(net.stats.melt, net.stats.rnd.melt)
# ggsave(filename=file.path("plots", "network stats - observed and simulated - reciprocity and assortment.png"), width=300, height=100, units="mm")


##  
flog.info("SNA: finished")
