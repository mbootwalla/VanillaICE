library(cacheSweave)
library(VanillaICE)
library(tools)
VanillaICE:::Sweave2pdf("CrlmmDownstream", driver=cacheSweaveDriver)
