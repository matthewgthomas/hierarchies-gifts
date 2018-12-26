##
## Analysis code for "Similar social complexity and cooperation in independent pastoralist societies"
##
## Author: Matthew Gwynfryn Thomas (@matthewgthomas)
##
## Packages used:
## - futile.logger
## - tidyverse
## - stringr
## - forcats
## - data.table
## - Matrix
## - RColorBrewer
## - rstanarm
## - loo (currently the development version for the model_weights() function)
## - bayesplot
## - xlsx
## - ineq
##
rm(list=ls())

# set up log file
library(futile.logger)
flog.appender(appender.file(paste0("analysis - ", format(Sys.time(), "%Y-%m-%d %H_%M_%S"), ".log")))  # create log file

# load all things
source("init.r")
source("init plots.r")
source("model comparison functions.r")

source("descriptive stats.r")

# social network analyses
source("SNA - reciprocity and assortment.r")
source("SNA - position and herd size.r")

# models predicting gift giving
source("multilevel models - fit.r")      # warning: this will take at least a week to fit all models (!)
source("multilevel models - compare.r")  # and this will take about 6 hours or so
source("multilevel models - analyse.r")

rmarkdown::render("results text.r", "word_document", output_dir=results.dir)  # produce the in-text statistics
rmarkdown::render("map.Rmd", output_dir=docs.dir, output_file="index.html")   # make the interactive map webpage

flog.info("Finished analysis!")
