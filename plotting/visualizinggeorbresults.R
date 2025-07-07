library(RevGadgets)
library(coda)
library(ggplot2)
library(ggtree)
library(grid)
library(gridExtra)

setwd("~/../../outputfilename")
file<-"simulationoutput.model.log"


trace_quant <- readTrace(path = file, burnin = 0.1)

summarizeTrace(trace = trace_quant, vars =  c("epoch_times[1]","epoch_times[2]"))

plotTrace(trace = trace_quant, vars = c("epoch_times[1]","epoch_times[2]"))[[1]]
