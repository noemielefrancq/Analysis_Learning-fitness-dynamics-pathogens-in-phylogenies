########################################################################################################################################
## Compute index dynamics pertussis
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

########################################################################################################################################
## Useful functions
########################################################################################################################################
## Load index functions
source('../Phylowave_Learning-fitness-dynamics-pathogens-in-phylogenies/2_Functions/2_1_Index_computation_20231220.R')

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
## Timed-trees, incorporating node support
## France
tree_France = read.nexus('3_BordetellaPertussis/1_Data/2_data/BEAST/hidden/France_mcc.tree')
data_beast = treeio::read.beast('3_BordetellaPertussis/1_Data/2_data/BEAST/hidden/France_mcc.tree')
data_beast = as.data.frame(data_beast@data)
data_beast = data_beast[which(as.numeric(data_beast$node) > length(tree_France$tip.label)),]
data_beast = data_beast[order(as.numeric(data_beast$node)),]
tree_France$node.label = as.numeric(data_beast$posterior)
tree_France = ladderize(multi2di(tree_France, random = F), right = F)

## SNP data
data_vcf_France = read.csv(file = "3_BordetellaPertussis/1_Data/2_data/FR_80p_dp5_SOR_QC_names_with_Tohama_recom_masked_snp_withN.vcf", sep = '\t')
colnames(data_vcf_France) = str_replace_all(colnames(data_vcf_France), '_', '\\|')
colnames(data_vcf_France)[which(colnames(data_vcf_France) == "NC|002929.2")] = "NC_002929.2"
pos_France = as.numeric(unlist(lapply(data_vcf_France$POS, function(x)str_split(x, '_')[[1]][1])))
for(i in 1:length(pos_France)){ ## Replace N and - by NAs in the vcf
  if(data_vcf_France$ALT[i] == 'N' |  data_vcf_France$ALT[i] == '-'){
    idx = which(pos_France == pos_France[i])
    idx = idx[which(data_vcf_France$ALT[idx] == 'A' | data_vcf_France$ALT[idx] == 'C'| 
                      data_vcf_France$ALT[idx] == 'T' | data_vcf_France$ALT[idx] == 'G')]
    for(j in 1:length(idx)){
      tmp = which(data_vcf_France[i,] == 1)
      data_vcf_France[idx[j],tmp] = NA
    }
  }
}
data_vcf_France = data_vcf_France[which(data_vcf_France$ALT != 'N'),]
data_vcf_France = data_vcf_France[which(data_vcf_France$ALT != '-'),]

## Metadata 
metadata_France = read.csv2('3_BordetellaPertussis/1_Data/2_data/France_strain_typing_202309.csv')

## Sort ptxP alleles: ptxp3 = 1, ptxp1 = 0, others = NA
metadata_France$ptxP[which(metadata_France$ptxP != 'ptxP1' & metadata_France$ptxP != 'ptxP3')] = NA
metadata_France$ptxP[which(metadata_France$ptxP == 'ptxP1')] = 0
metadata_France$ptxP[which(metadata_France$ptxP == 'ptxP3')] = 1

## Sort fim3 alleles: fim3-1 = 1, fim3-2 = 1,others = NA
metadata_France$fim3[which(metadata_France$fim3 != 'fim3-1' & metadata_France$fim3 != 'fim3-2')] = NA
metadata_France$fim3[which(metadata_France$fim3 == 'fim3-1')] = 0
metadata_France$fim3[which(metadata_France$fim3 == 'fim3-2')] = 1

## Names all sequences
names_seqs_France = tree_France$tip.label

## Mutation rate from BEAST
mu = 2.5E-7 ## in mutations/site/year

## Length genome 
genome_length = 4086189 ## Length genome

## Parameters for THD
timescale = 2 ## Timescale
bandwidth = thd.bandwidth(timescale, genome_length, mu, q = 0.5) ## Corresponding bandwidth
########################################################################################################################################

########################################################################################################################################
## Preparation data tips
########################################################################################################################################
## Create dataset with names of each sequence, time sampling, prn type etc
dataset_tips_France = data.frame('ID' = 1:length(names_seqs_France),
                                 'name_seq' = names_seqs_France,
                                 'time' = unlist(lapply(names_seqs_France, function(x)str_split(x, pattern = "\\|")[[1]][2])))
dataset_tips_France$time[which(dataset_tips_France$name_seq == 'NC_002929.2')] = 1954
dataset_tips_France$time = as.numeric(dataset_tips_France$time)
## Add SNP data to the main dataset
a = match(names_seqs_France, colnames(data_vcf_France))
dataset_tips_France = cbind(dataset_tips_France, t(data_vcf_France[,a]))
colnames(dataset_tips_France) = c('ID', 'name_seq', 'time', data_vcf_France$POS)
########################################################################################################################################

########################################################################################################################################
## Preparation data nodes
########################################################################################################################################
## Compute distance between each pair of sequences AND NODES in the tree
genetic_distance_mat_France = dist.nodes.with.names(tree_France)

## Get the time each node
nroot = length(tree_France$tip.label)+1 ## Checked and it's the root 
distance_to_root = genetic_distance_mat_France[nroot,]
root_height_France = dataset_tips_France$time[which(dataset_tips_France$name_seq == names(distance_to_root[20]))] - distance_to_root[20]  ## Take one tip, doesn't matter which tip is used
nodes_height = root_height_France + distance_to_root[length(names_seqs_France)+(1:(length(names_seqs_France)-1))]

# Meta-data with nodes 
dataset_with_nodes_France = data.frame('ID' = c(1:length(names_seqs_France), length(names_seqs_France)+(1:(length(names_seqs_France)-1))),
                                       'name_seq' = c(names_seqs_France, length(names_seqs_France)+(1:(length(names_seqs_France)-1))),
                                       'time' = c(dataset_tips_France$time, nodes_height),
                                       'is.node' = c(rep('no', length(names_seqs_France)), rep('yes', (length(names_seqs_France)-1)))) 
m = match(dataset_with_nodes_France$name_seq, metadata_France$SEQ)
idx = which(is.na(m) == F)
m = m[idx]
dataset_with_nodes_France$genotype = NA
dataset_with_nodes_France$genotype[idx] = paste0(metadata_France$ptxP[m], '_', metadata_France$fim3[m])
dataset_with_nodes_France$genotype[idx[which(metadata_France$ptxP[m] == 0 & is.na(metadata_France$fim3[m]))]] = 'NA_NA' 
dataset_with_nodes_France$genotype[idx[which(metadata_France$ptxP[m] != 0 & is.na(metadata_France$fim3[m]))]] = 'NA_NA'
dataset_with_nodes_France$genotype[idx[which(is.na(metadata_France$ptxP[m]))]] = 'NA_NA'
dataset_with_nodes_France$genotype[which(dataset_with_nodes_France$name_seq == "NC_002929.2")] = '0_0'
########################################################################################################################################

########################################################################################################################################
## Reconstruction genotypes at each node
########################################################################################################################################
## France
snp_data = dataset_with_nodes_France$genotype[which(dataset_with_nodes_France$is.node == 'no')]
snp_data = as.numeric(as.factor(snp_data))
snp_data[which(snp_data == 5)] = NA
tree = tree_France
tree$node.label = NULL
rec = ace(snp_data, tree, type = 'discrete', CI = T)
rec_nodes = apply(rec$lik.anc, MAR = 1, function(x)which.max(x))
rec_all = c(as.numeric(as.character(snp_data)), rec_nodes) ## List of all states: first all tips, then all nodes
dataset_with_nodes_France$genotype = rec_all
########################################################################################################################################

########################################################################################################################################
## Compute index of every tip and node
########################################################################################################################################
## Window of time on which to search for samples in the population
wind = 1 ## years

dataset_with_nodes_France$index = compute.index(time_distance_mat = genetic_distance_mat_France, 
                                                       timed_tree = tree_France, 
                                                       time_window = wind,
                                                       metadata = dataset_with_nodes_France, 
                                                       mutation_rate = mu,
                                                       timescale = timescale,
                                                       genome_length = genome_length)
########################################################################################################################################

########################################################################################################################################
## Generate colors needed
########################################################################################################################################
library(MetBrewer)
lev = as.character(1:4)
n = length(lev)
colors_genotypes = c(met.brewer(name="Cross", n=7, type="continuous")[c(1,5,6)], 'grey60')
colors_genotypes_2 = c('#cccccc', '#666666', '#333333', 'white')

dataset_with_nodes_France$genotype_color = dataset_with_nodes_France$genotype
dataset_with_nodes_France$genotype_color = as.factor(dataset_with_nodes_France$genotype_color)
labels = levels(dataset_with_nodes_France$genotype_color)
levels(dataset_with_nodes_France$genotype_color) = colors_genotypes[match(labels, lev)]
dataset_with_nodes_France$genotype_color = as.character(dataset_with_nodes_France$genotype_color)
dataset_with_nodes_France$genotype_color[which(is.na(dataset_with_nodes_France$genotype_color))] = 'grey60'
########################################################################################################################################

########################################################################################################################################
## Save output
########################################################################################################################################
save.image(file='3_BordetellaPertussis/2_analysis_index/1_index_computations/Initial_index_computation_and_parameters.Rdata')
########################################################################################################################################

