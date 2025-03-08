########################################################################################################################################
## Detect lineages for each country, SC2
########################################################################################################################################
## Packages used
library(stringr)
library(ape)
library(phangorn)
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

########################################################################################################################################
## Load data
########################################################################################################################################
load('1_SARS-CoV-2/2_analysis_index/1_index_computations/Initial_index_computation_and_parameters.Rdata')

########################################################################################################################################
## Useful functions
########################################################################################################################################
## Load index functions
source('../Phylowave_Learning-fitness-dynamics-pathogens-in-phylogenies/2_Functions/2_1_Index_computation_20231220.R')
source('../Phylowave_Learning-fitness-dynamics-pathogens-in-phylogenies/2_Functions/2_2_Lineage_detection_20231220.R')

## Other functions
mean.and.ci <-function(v){ return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))}
axisPhylo_NL = function (side = 1, root.time = NULL, backward = TRUE, at_axis = NULL, lab_axis = NULL, ...){
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  type <- lastPP$type
  if (type == "unrooted")
    stop("axisPhylo() not available for unrooted plots; try add.scale.bar()")
  if (type == "radial")
    stop("axisPhylo() not meaningful for this type of plot")
  if (is.null(root.time))
    root.time <- lastPP$root.time
  if (type %in% c("phylogram", "cladogram")) {
    xscale <- if (lastPP$direction %in% c("rightwards", "leftwards"))
      range(lastPP$xx)
    else range(lastPP$yy)
    tmp <- lastPP$direction %in% c("leftwards", "downwards")
    tscale <- c(0, xscale[2] - xscale[1])
    if (xor(backward, tmp))
      tscale <- tscale[2:1]
    if (!is.null(root.time)) {
      tscale <- tscale + root.time
      if (backward)
        tscale <- tscale - xscale[2]
    }
    beta <- diff(xscale)/diff(tscale)
    alpha <- xscale[1] - beta * tscale[1]
    if(is.null(at_axis) == T){
      x <- beta * lab + alpha
      lab <- pretty(tscale)
    }
    # if(is.null(at_axis) != F){
    x <- at_axis
    lab <- lab_axis
    # }
    axis(side = side, at = x, labels = lab, ...)
  }
  else {
    n <- lastPP$Ntip
    xx <- lastPP$xx[1:n]
    yy <- lastPP$yy[1:n]
    r0 <- max(sqrt(xx^2 + yy^2))
    alpha <- sort(setNames(rect2polar(xx, yy)$angle, 1:n))
    angles <- c(diff(alpha), 2 * pi - alpha[n] + alpha[1L])
    j <- which.max(angles)
    i <- if (j == 1L)
      n
    else j - 1L
    firstandlast <- as.integer(names(angles[c(i, j)]))
    theta0 <- mean(atan2(yy[firstandlast], xx[firstandlast]))
    x0 <- r0 * cos(theta0)
    y0 <- r0 * sin(theta0)
    inc <- diff(pretty(c(0, r0))[1:2])
    srt <- 360 * theta0/(2 * pi)
    coef <- -1
    if (abs(srt) > 90) {
      srt <- srt + 180
      coef <- 1
    }
    len <- 0.025 * r0
    r <- r0
    while (r > 1e-08) {
      x <- r * cos(theta0)
      y <- r * sin(theta0)
      if (len/r < 1) {
        ra <- sqrt(len^2 + r^2)
        thetaa <- theta0 + coef * asin(len/r)
        xa <- ra * cos(thetaa)
        ya <- ra * sin(thetaa)
        segments(xa, ya, x, y)
        text(xa, ya, r0 - r, srt = srt, adj = c(0.5,
                                                1.1), ...)
      }
      r <- r - inc
    }
    segments(x, y, x0, y0)
  }
}
read.chains.from.table = function(table){
  if(typeof(table) != 'double') {
    print('Changing type of table to matrix')
    table = as.matrix(table)
  }
  Chains = list()
  col_variables = sapply(colnames(table), function(x)str_split(x, pattern = "[.]")[[1]][1])
  variable_names = unique(col_variables)
  nchains = nrow(table)
  for(i in 1:length(variable_names)){
    a = match(col_variables, variable_names[i])
    if(length(which(is.na(a) == F)) == 1) {
      Chains[[i]] = table[,match(variable_names[i], col_variables)]
    }
    else{
      a = match(col_variables, variable_names[i])
      tmp = colnames(table)[which(is.na(a) == F)]
      ndims = 0
      dims = NULL
      empty = F
      while(empty == F & ndims < 10){
        l = sapply(tmp, function(x)str_split(x, pattern = "[.]")[[1]][1+ndims+1])
        if(all(is.na(l))) empty = T
        if(all(is.na(l)) == F) {
          ndims = ndims + 1
          dims = c(dims, max(as.numeric(l)))
        }
      }
      if(ndims >5) print('Error: this function only supports arrays of <=5 dimensions')
      Chains[[i]] = array(NA, dim = c(nchains, dims))
      a = which(is.na(match(col_variables, variable_names[i]))==F)
      
      if(ndims == 1){
        Chains[[i]] = table[,a]
        colnames(Chains[[i]]) = NULL
      }
      
      if(ndims == 2){
        j = 1
        for(d2 in 1:dims[2]){
          for(d1 in 1:dims[1]){
            Chains[[i]][,d1,d2] = table[,a[j]]
            j = j+1
          }
        }
      }
      
      if(ndims == 3){
        j = 1
        for(d3 in 1:dims[3]){
          for(d2 in 1:dims[2]){
            for(d1 in 1:dims[1]){
              Chains[[i]][,d1,d2,d3] = table[,a[j]]
              j = j+1
            }
          }
        }
      }
      
      if(ndims == 4){
        j = 1
        for(d4 in 1:dims[4]){
          for(d3 in 1:dims[3]){
            for(d2 in 1:dims[2]){
              for(d1 in 1:dims[1]){
                Chains[[i]][,d1,d2,d3,d4] = table[,a[j]]
                j = j+1
              }
            }
          }
        }
      }
      
      if(ndims == 5){
        j = 1
        for(d5 in 1:dims[5]){
          for(d4 in 1:dims[4]){
            for(d3 in 1:dims[3]){
              for(d2 in 1:dims[2]){
                for(d1 in 1:dims[1]){
                  Chains[[i]][,d1,d2,d3,d4] = table[,a[j]]
                  j = j+1
                }
              }
            }
          }
        }
      }
    }
    names(Chains)[i] = variable_names[i]
  }
  return(Chains)
}
########################################################################################################################################

######################################################################################################################################
## Read best splits
######################################################################################################################################
k_smooth = -1
min_group_size = 2
potential_splits_world = readRDS(paste0('1_SARS-CoV-2/2_analysis_index/2_find_index_groups/potential_splits_weighting_world_tracked_ksmooth', k_smooth, '_min_group_size_', min_group_size, '.rds'))$best_nodes_names[[15]]
split_world = merge.groups(timed_tree = tree_sars_world, metadata = dataset_with_nodes_world, 
                            initial_splits = potential_splits_world, 
                            group_count_threshold = 10, group_freq_threshold = 0)
######################################################################################################################################

#######################################################################################################################################
## Assign colour to each group
#######################################################################################################################################
library(MetBrewer)
## World
dataset_with_nodes_world$groups = as.factor(split_world$groups)
## Reorder labels by time of emergence
name_groups_world = levels(dataset_with_nodes_world$groups)
time_groups_world = NULL
for(i in 1:length(name_groups_world)){
  time_groups_world = c(time_groups_world, min(dataset_with_nodes_world$time[which(dataset_with_nodes_world$groups == name_groups_world[i] &
                                                                                   dataset_with_nodes_world$is.node == 'no')]))
}
levels(dataset_with_nodes_world$groups) = match(name_groups_world, order(time_groups_world, decreasing = T))
dataset_with_nodes_world$groups = as.numeric(as.character(dataset_with_nodes_world$groups))
dataset_with_nodes_world$groups = as.factor(dataset_with_nodes_world$groups)
## Update split
split_world$tip_and_nodes_groups = match(split_world$tip_and_nodes_groups, order(time_groups_world, decreasing = T))
names(split_world$tip_and_nodes_groups) = 1:length(split_world$tip_and_nodes_groups)
split_world$groups = as.factor(split_world$groups)
levels(split_world$groups) = match(name_groups_world, order(time_groups_world, decreasing = T))
split_world$groups = as.numeric(as.character(split_world$groups))
## Choose color palette
n_groups_world <- length(name_groups_world)
colors_groups_world = (met.brewer(name="Cross", n=n_groups_world, type="continuous"))
## Color each group
dataset_with_nodes_world$group_color = dataset_with_nodes_world$groups
levels(dataset_with_nodes_world$group_color) = colors_groups_world
dataset_with_nodes_world$group_color = as.character(dataset_with_nodes_world$group_color)
#######################################################################################################################################

#######################################################################################################################################
## Save the detection
#######################################################################################################################################
save.image('1_SARS-CoV-2/2_analysis_index/2_find_index_groups/Lineages_detected.Rdata')
#######################################################################################################################################

