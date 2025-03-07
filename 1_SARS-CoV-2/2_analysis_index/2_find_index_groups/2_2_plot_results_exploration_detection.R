#############################################################
## Plot results detection
#############################################################
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
library(cowplot)
library(scales)
font_import()
loadfonts(device="all")

#############################################################
## Read results
#############################################################
k_smooth = -1
min_group_size = 2

potential_splits_world = readRDS(paste0('1_SARS-CoV-2/2_analysis_index/2_find_index_groups/potential_splits_weighting_world_tracked_ksmooth', k_smooth, '_min_group_size_', min_group_size, '.rds'))
#############################################################

#############################################################
## Format in results in dataframes
#############################################################
world_explained_dev = data.frame('N_groups' = 0:length(potential_splits_world$best_dev_explained),
                                     'Non_explained_deviance' = (1-c(potential_splits_world$first_dev, potential_splits_world$best_dev_explained)),
                                     'Non_explained_deviance_log' = log(1-c(potential_splits_world$first_dev, potential_splits_world$best_dev_explained)))
world_explained_dev$Non_explained_deviance_log = world_explained_dev$Non_explained_deviance_log-min(world_explained_dev$Non_explained_deviance_log)
world_explained_dev = world_explained_dev[1:31,]
#############################################################

#############################################################
## world
#############################################################
pdf('1_SARS-CoV-2/2_analysis_index/2_find_index_groups/Elbow_plot_explained_deviance_N_groups_wolrd.pdf', 
    width = 10/2.54, height = 5/2.54)
par(mfrow = c(1,2), oma = c(0.5,0.5,0.5,0.5), mar = c(1.6,1.6,0,0.5), mgp = c(0.7,0,0), family = 'Arial',
    cex.axis=0.5, cex.lab=0.5, cex.main=0.7, cex.sub=0.5)
# par(mfrow = c(1,1))
plot(world_explained_dev$N_groups,
     world_explained_dev$Non_explained_deviance,
     bty = 'n', ylim = c(0, ceiling(10*max(world_explained_dev$Non_explained_deviance))/10),
     yaxt = 'n', pch = 16, main = 'world - Linear', cex = 0.5,
     xaxt = 'n', ylab = 'Non-explained deviance (%)', xlab = 'Number of groups')
axis(1, lwd = 0.5, tck=-0.02)
axis(2, las = 2, at = seq(0,ceiling(10*max(world_explained_dev$Non_explained_deviance))/10,0.1),
     labels = seq(0,ceiling(10*max(world_explained_dev$Non_explained_deviance))/10,0.1)*100, lwd = 0.5, tck=-0.02)
abline(v = 15, lwd = 2, lty = 2, col = 'chartreuse4')

plot(world_explained_dev$N_groups,
     (world_explained_dev$Non_explained_deviance),
     ylim = c(0.001, 1),
     bty = 'n', log = 'y', cex = 0.5,
     yaxt = 'n', pch = 16, main = 'world - log',
     xaxt = 'n', ylab = 'Non-explained deviance (%)', xlab = 'Number of groups')
axis(1, lwd = 0.5, tck=-0.02)
axis(2, las = 2, at = 10^seq(-3,0,1),
     labels = 10^seq(-3,0,1)*100, lwd = 0.5, tck=-0.02)
abline(v = 15, lwd = 2, lty = 2, col = 'chartreuse4')

## Chosen iteration
## 15
## dev expalined = 0.9933211
dev.off()
#############################################################

