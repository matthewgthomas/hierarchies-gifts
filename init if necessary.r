##
## A kind of meta-loader to import data if not already in the workspace
## This lets you run each part of the analysis separately rather than all in a batch, should you prefer
##
if (!exists("dyads"))
  source("init.r")

if (!exists("common_theme"))
  source("init plots.r")
