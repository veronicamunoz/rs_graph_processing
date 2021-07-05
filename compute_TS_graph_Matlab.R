#!/usr/bin/env Rscript
library('bitops')
library('tools')
library('oro.nifti')
library('R.matlab')
print("OK")
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  
  stop("ERROR. No arguments supplied.n", call.=FALSE)
  
} else if (length(args) == 7) { # Compute only time series 
  setwd(args[1])
  source('computeTS_annotated.R')
  computeTS_annotated(rs_path = args[2], 
                      r_dyn_name = args[3], 
                      art_path = args[4], 
                      atlas_file = args[5], 
                      gm_file = args[6], 
                      graph_path = args[7]
  )
} else if (length(args) == 11) { # Compute only graph analysis (time series have already been calculated)
  print("hello")
  setwd(args[1])
  source('Figure.R')
  source('computeGraphs_annotated.R')
  compute_Graph_annotated(rs_path = args[2],
                          atlas_file = args[3],
                          graph_path = args[4],
                          file_coord = args[5],
                          graphs = args[6],
                          regions_selected = eval(parse(text = args[7])),
                          percentage_selected = eval(parse(text = args[8])),
                          num.levels = eval(parse(text = args[9])),
                          num.nb.edges = eval(parse(text = args[10])),
                          length_time_series = eval(parse(text = args[11]))
  )
  
} 





  