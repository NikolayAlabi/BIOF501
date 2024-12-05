#! /usr/bin/env Rscript
library(ape)

args<-commandArgs(TRUE)

filename <- args[1]
plotname <- args[2]

tree <- ape::read.tree(filename)

png(plotname)
plot(tree)
dev.off
