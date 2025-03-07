########################################################################################################################################
## Comparison detected groups vs. genotypes
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
library(cowplot)
library(extrafont)
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
## Load data
########################################################################################################################################
load('1_SARS-CoV-2/2_analysis_index/2_find_index_groups/Lineages_detected.Rdata')
########################################################################################################################################

########################################################################################################################################
## Create VOC column
########################################################################################################################################
dataset_with_nodes_world$VOC = dataset_with_nodes_world$clade
dataset_with_nodes_world$VOC = factor(dataset_with_nodes_world$VOC)
levels(dataset_with_nodes_world$VOC) = c('wt', 'wt', 'wt', 'wt', 'wt', 'wt', "20E/B.1.177 - EU1", 'wt', 'wt', "20H/B.1.351 - Beta", "20I/B.1.1.7 - Alpha", 
                                         "20J/P.1 - Gamma", "21A-I-J /B.1.617.2 - Delta", "21B/B.1.617.1 - Kappa", "21C/B.1.427 - Epsilon", "21D/B.1.525 - Eta", "21F/B.1.526 - Iota", "21G/C.37 - Lambda",
                                         "21H/B.1.621 - Mu", "21A-I-J /B.1.617.2 - Delta", "21A-I-J /B.1.617.2 - Delta", "21K / BA.1 - Omicron", "21L / ~BA.2 - Omicron","22A / BA.4 - Omicron", "22B / BA.5 - Omicron", 
                                         "22C / BA.2.12.1 - Omicron", "22D / BA.2.75 - Omicron", "22E / BQ.1 - Omicron", "22F / XBB - Omicron", "23A / XBB.1.5 - Omicron")
########################################################################################################################################

########################################################################################################################################
## Find colors Nextstrain clades and VOC
########################################################################################################################################
colors_groups_world_nextstrain = MetBrewer::met.brewer(name="Hiroshige", n=length(levels(as.factor(dataset_with_nodes_world$clade))), type="continuous")
names(colors_groups_world_nextstrain) = levels(as.factor(dataset_with_nodes_world$clade))
colors_groups_world_VOC = c(MetBrewer::met.brewer(name="Hiroshige", n=length(levels(as.factor(dataset_with_nodes_world$VOC))), type="continuous"))
names(colors_groups_world_VOC) = levels(as.factor(dataset_with_nodes_world$VOC))
########################################################################################################################################

#######################################################################################################################################
## Compare groups with VOCs
#######################################################################################################################################
groups = matrix(dataset_with_nodes_world$groups[which(dataset_with_nodes_world$is.node == 'no')], ncol = 1)
colnames(groups) = 'groups'
rownames(groups) = dataset_with_nodes_world$name_seq[which(dataset_with_nodes_world$is.node == 'no')]
cols = colors_groups_world
names(cols) = as.character(1:max(as.numeric(name_groups_world)))
plot_tree_sars_world_groups <- ggtree(tree_sars_world, mrsd=lubridate::date_decimal(max(times_seqs_world)), size = 0.08,
                                      aes(color = as.character(dataset_with_nodes_world$groups))) + 
  scale_color_manual(values = cols)+
  theme_tree2()
plot_tree_sars_world_groups = gheatmap(plot_tree_sars_world_groups, groups, offset=0.1, width=0.10, 
                                       colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = cols)+
  # scale_x_reverse() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

VOC = matrix(dataset_with_nodes_world$VOC[which(dataset_with_nodes_world$is.node == 'no')], ncol = 1)
colnames(VOC) = 'groups'
rownames(VOC) = dataset_with_nodes_world$name_seq[which(dataset_with_nodes_world$is.node == 'no')]
plot_tree_sars_world_VOC <- ggtree(tree_sars_world, mrsd=lubridate::date_decimal(max(times_seqs_world)),size = 0.08,
                                   aes(color = as.character(dataset_with_nodes_world$VOC))) + 
  scale_color_manual(values = colors_groups_world_VOC)+
  theme_tree2(legend = 'none')
plot_tree_sars_world_VOC = gheatmap(plot_tree_sars_world_VOC, VOC, offset=0.1, width=0.10, 
                                           colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = colors_groups_world_VOC, na.value = 'white')+
  scale_x_reverse() +
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

p = plot_grid(plot_tree_sars_world_groups, plot_tree_sars_world_VOC, 
          rel_widths = c(1, 1), labels = '', ncol = 2)
ggsave(filename = '1_SARS-CoV-2/2_analysis_index/2_find_index_groups/Comparison_groups_VOC_world.pdf',
       plot = p, device = 'pdf', scale = 1, width = 7, height = 7, units = 'cm')
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs known lineages (world, scale)
#######################################################################################################################################
## Remove clades that are present in only tiny numbers
dataset_with_nodes_world$VOC = as.character(dataset_with_nodes_world$VOC)
tmp = names(which(table(dataset_with_nodes_world$VOC) < 0))
dataset_with_nodes_world$VOC[which(!is.na(match(dataset_with_nodes_world$VOC, tmp)))] = NA
# Reset levels
dataset_with_nodes_world$VOC = as.factor(as.character(dataset_with_nodes_world$VOC))

## Create VOC column removing small clades
grouping1 = levels(as.factor(dataset_with_nodes_world$VOC))
grouping2 = levels(as.factor(dataset_with_nodes_world$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes_world$VOC == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes_world$groups[idx])))
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = TRUE,
  squared_order = T,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}
#######################################################################################################################################

#######################################################################################################################################
## ARI world
#######################################################################################################################################
library(fossil)
idx = which(!is.na(dataset_with_nodes_world$VOC))
group1 = as.numeric(dataset_with_nodes_world$VOC[idx])
group2 = dataset_with_nodes_world$groups[idx]
rand.treeannotator_partitions = fossil::adj.rand.index(group1, group2)
rand.treeannotator_partitions
# 0.8021831
#######################################################################################################################################
