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
library(cowplot)
library(scales)
font_import()
loadfonts(device="all")

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

########################################################################################################################################
## Load data & set parameters
########################################################################################################################################
## Timed-tree
tree_sars_world = collapse.singles(ladderize(multi2di(read.nexus('1_SARS-CoV-2/1_Data_Nextstrain_20230414/World/nextstrain_ncov_gisaid_global_all-time_timetree_tipdates_tipclades.nexus'), random = F)))

## Names all sequences
names_seqs_world = tree_sars_world$tip.label

## Collection times of all sequences
times_seqs_world = as.numeric(sapply(names_seqs_world, function(x)tail(str_split(x, pattern = '/')[[1]],2)[1]))

## Nextstrain clades of all sequences
clades_seqs_world = sapply(names_seqs_world, function(x)tail(str_split(x, pattern = '/')[[1]],1))

## Index parameters
## Length genome 
genome_length = 29903 # reference nextstrain https://www.ncbi.nlm.nih.gov/nuccore/MN908947
## Mutation rate 
mutation_rate = 8.1e-4 # mutation rate used by nextstrain https://github.com/nextstrain/ncov

## Parameters for the index
timescale = 0.15 ## Timescale
# timescale = 20 ## Timescale

## Window of time on which to search for samples in the population
wind = 15 #days
wind = wind/365
########################################################################################################################################

########################################################################################################################################
## Genetic distance
########################################################################################################################################
genetic_distance_mat_world = dist.nodes.with.names(tree_sars_world)
########################################################################################################################################

########################################################################################################################################
## Get the time of each node & prepare data tips and nodes
########################################################################################################################################
## World
nroot= length(tree_sars_world$tip.label) + 1 ## Root number
distance_to_root = genetic_distance_mat_world[nroot,]
root_height_world = times_seqs_world[which(names_seqs_world == names(distance_to_root[1]))] - distance_to_root[1]
nodes_height = root_height_world + distance_to_root[length(names_seqs_world)+(1:(length(names_seqs_world)-1))]
dataset_with_nodes_world = data.frame('ID' = c(1:length(names_seqs_world), length(names_seqs_world)+(1:(length(names_seqs_world)-1))),
                                      'name_seq' = c(names_seqs_world, length(names_seqs_world)+(1:(length(names_seqs_world)-1))),
                                      'time' = c(times_seqs_world, nodes_height),
                                      'is.node' = c(rep('no', length(names_seqs_world)), rep('yes', (length(names_seqs_world)-1))),
                                      'clade' = c(clades_seqs_world, rep(NA, length(names_seqs_world)-1)))
########################################################################################################################################

########################################################################################################################################
## Reconstruction clade at each node
########################################################################################################################################
## World
snp_data = dataset_with_nodes_world$clade[which(dataset_with_nodes_world$is.node == 'no')]
snp_data = as.factor(snp_data)
lev = levels(snp_data)
tree = tree_sars_world
tree$node.label = NULL
rec = ace(snp_data, tree, type = 'discrete', CI = T)
rec_nodes = apply(rec$lik.anc, MAR = 1, function(x)which.max(x))
rec_all = c(as.numeric(snp_data), rec_nodes) ## List of all states: first all tips, then all nodes
rec_all = as.factor(rec_all)
levels(rec_all) = lev
dataset_with_nodes_world$clade = as.character(rec_all)
########################################################################################################################################

########################################################################################################################################
## Compute index
########################################################################################################################################
dataset_with_nodes_world$index = compute.index(time_distance_mat = genetic_distance_mat_world, 
                                               timed_tree = tree_sars_world, 
                                               time_window = wind,
                                               metadata = dataset_with_nodes_world, 
                                               mutation_rate = mutation_rate,
                                               timescale = timescale,
                                               genome_length = genome_length)
########################################################################################################################################

########################################################################################################################################
## Colors datasets, given nextstrain clades
########################################################################################################################################
library(MetBrewer)
all_clades = c(dataset_with_nodes_world$clade)
clades = c(levels(as.factor(all_clades)), NA)
colors_clade = met.brewer(name="Cross", n=length(levels(as.factor(all_clades))), type="continuous")

## World
dataset_with_nodes_world$color = factor(dataset_with_nodes_world$clade, levels = clades)
levels(dataset_with_nodes_world$color) = colors_clade
dataset_with_nodes_world$color = as.character(dataset_with_nodes_world$color)
dataset_with_nodes_world$color[which(is.na(dataset_with_nodes_world$color))] = 'grey'
########################################################################################################################################

########################################################################################################################################
## Save output
########################################################################################################################################
save.image(file='1_SARS-CoV-2/2_analysis_index/1_index_computations/Initial_index_computation_and_parameters.Rdata')
########################################################################################################################################







