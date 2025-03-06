########################################################################################################################################
## Initial plots index Pertussis, by genotype
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

load('3_BordetellaPertussis/2_analysis_index/1_index_computations/Initial_index_computation_and_parameters.Rdata')

########################################################################################################################################
## Plot France
########################################################################################################################################
pdf('3_BordetellaPertussis/2_analysis_index/1_index_computations/Index_pertussis_france_genotypes.pdf', width = 5/2.54, height = 5/2.54)
par(mfcol = c(2,1), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)

min_year = 1970
max_year = 2023
max_index = 0.45

## France
tree = tree_France
root_height = root_height_France
dataset_with_nodes = dataset_with_nodes_France
## Plot tree
tree_to_plot = tree
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.2, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$genotype_color[which(dataset_with_nodes$is.node == 'no')], cex = 0.2)
# nodelabels(pch = 16, col = dataset_with_nodes$genotype_color[which(dataset_with_nodes$is.node == 'yes')], cex = 0.5)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 10)-root_height,
             lab_axis = seq(min_year, max_year, 10), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')], 
     dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes$genotype_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.2, 
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
       dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes$genotype_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.2, pch = 16)
axis(1, at = seq(min_year, max_year, 10), labels = seq(min_year, max_year, 10), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="Pertussis - France", line=-0.5, outer = F)
title(ylab="Diversity index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)
dev.off()
########################################################################################################################################

