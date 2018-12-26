##
## Options, parameters, etc. common to all plots
##
library(ggplot2)

cols = c("#5F8A71", "#2D4A94")  # blue and green

# saami flag colours
saami_blue_green = c("#5F8A71", "#2D4A94")
saami_red = rgb(193, 48, 72, 255, maxColorValue=255)
saami_lightred = rgb(206, 99, 122,  255, maxColorValue=255)
saami_yellow = rgb(227, 204, 134, maxColorValue=255)
saami_yellow_dark = rgb(205, 163, 37, maxColorValue = 255)


common_theme = theme_bw() +
  #eliminates baground, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
    ,axis.text=element_text(size=12) 
    ,axis.title=element_text(size=12)  # 12 for printing; 22 for poster
  ) +
  theme(legend.position="none")

poster_theme = theme(
     axis.text=element_text(size=28)  # size=12 for printing; 28 for poster
    ,axis.title=element_text(size=28)  # 12 for printing; 28 for poster
    ,strip.text.x = element_text(size=24)  # facet text size
  )
