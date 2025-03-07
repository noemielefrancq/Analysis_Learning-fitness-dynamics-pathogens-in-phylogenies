########################################################################################################################################
## Initial plots index SC2, by nextrsain clade
########################################################################################################################################
## Packages used
library(stringr)
library(ape)
library(phangorn)
library(BactDating)
library(phytools)
library(coda)
library(thd)
library(vcfR)
library(lubridate)
library(ggplot2)
library(ggtree)
library(extrafont)
# font_import()
loadfonts(device="all")

load('1_SARS-CoV-2/2_analysis_index/1_index_computations/Initial_index_computation_and_parameters.Rdata')

########################################################################################################################################
## Plot tree & THD below, with colors from genotypes
## WORLD
########################################################################################################################################
pdf('1_SARS-CoV-2/2_analysis_index/1_index_computations/Index_sarscov2_world_nextstrainclades.pdf', width = 5/2.54, height = 5/2.54)
par(mfcol = c(2,1), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)

min_year = 2020
max_year = 2023.5
max_index = 1

## World
tree = tree_sars_world
root_height = root_height_world
dataset_with_nodes = dataset_with_nodes_world
## Plot tree
tree_to_plot = ladderize(tree, right = F)
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$color, cex = 0.15)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 0.5)-root_height,
             lab_axis = seq(min_year, max_year, 0.5), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')], 
     dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes$color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.15, 
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
       dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes$color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.15, pch = 16)
axis(1, at = seq(min_year, max_year, 0.5), labels = seq(min_year, max_year, 0.5), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="SARS-CoV-2", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

## If leggend is needed
# idx = which(!duplicated(dataset_with_nodes_world$clade))
# legend('topleft', legend = dataset_with_nodes_world$clade[idx], col = dataset_with_nodes_world$color[idx], pch = 16, cex = 0.2)
dev.off()
########################################################################################################################################




