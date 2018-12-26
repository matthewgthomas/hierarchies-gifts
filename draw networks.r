##
## Draw networks
##
library(tidyverse)
library(ggraph)
library(igraph)
source("init if necessary.r")

for (site in unique(people$Pop)) {
  gifts = dyads %>% 
    filter(Pop == site, GiftGiven > 0) %>% 
    mutate(Weight = as.integer(as.factor(Amount)) / 10) %>% 
    select(Ego = Ego.ID, Alter = Alter.ID, Weight)
  
  nodes = people %>% 
    filter(Pop == site)
  
  g = graph_from_data_frame(gifts, vertices = nodes)
  
  gg = ggraph(g) +
    geom_edge_link(aes(edge_width = Weight),
                   colour = "lightgrey",
                   arrow = arrow(length = unit(1.5, 'mm')), 
                   end_cap = circle(3, 'mm')) + 
    geom_node_point(aes(size = gift.deg.in, colour = HG)) +
    scale_edge_width_continuous(range = c(0.5, 2)) +
    theme_void() +
    theme(legend.position="none")
  
  ggsave(paste0("plots/network - ", site, ".pdf"), plot = gg,
         width = 200, height = 200, units = "mm")
  
  print(paste0("Plotted ", site))
}
