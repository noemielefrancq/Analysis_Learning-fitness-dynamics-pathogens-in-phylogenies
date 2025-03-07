########################################################################################################################################
## Associate SNPs to groups
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
library(cowplot)
library(scales)
# font_import()
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
load('1_SARS-CoV-2/2_analysis_index/1_index_computations/Initial_index_computation_and_parameters.Rdata')
load('1_SARS-CoV-2/2_analysis_index/2_find_index_groups/Lineages_detected.Rdata')
########################################################################################################################################

########################################################################################################################################
## Load AA data
########################################################################################################################################
## Vcfs
########################################################################################################################################
ORFs = c('E', 'M', 'N', 'ORF10', 'ORF14', 
         'ORF1a', 'ORF1b', 'ORF3a', 'ORF6', 'ORF7a', 'ORF7b',
         'ORF8', 'ORF9b', 'S')

data_vcf_world_E = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_E.vcf", sep = '\t')
data_vcf_world_M = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_M.vcf", sep = '\t')
data_vcf_world_N = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_N.vcf", sep = '\t')
data_vcf_world_ORF10 = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF10.vcf", sep = '\t')
data_vcf_world_ORF14 = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF14.vcf", sep = '\t')
data_vcf_world_ORF1a = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF1a.vcf", sep = '\t')
data_vcf_world_ORF1b = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF1b.vcf", sep = '\t')
data_vcf_world_ORF3a = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF3a.vcf", sep = '\t')
data_vcf_world_ORF6 = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF6.vcf", sep = '\t')
data_vcf_world_ORF7a = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF7a.vcf", sep = '\t')
data_vcf_world_ORF7b = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF7b.vcf", sep = '\t')
data_vcf_world_ORF8 = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF8.vcf", sep = '\t')
data_vcf_world_ORF9b = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF9b.vcf", sep = '\t')
data_vcf_world_S = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_S.vcf", sep = '\t')
data_vcf_world_full = read.csv(file = "1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_full.vcf", sep = '\t')

## Reconstructions
idx_min = which.min(dataset_with_nodes_world$time[which(dataset_with_nodes_world$is.node == 'no')])

reconstruction_world_E = readRDS('1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_E.rds')
reconstruction_world_E$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_E$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_E$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_E$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_E$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_E$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_E$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_E$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_E$possible_snps = colnames(reconstruction_world_E$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_M = readRDS('1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_M.rds')
reconstruction_world_M$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_M$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_M$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_M$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_M$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_M$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_M$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_M$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_M$possible_snps = colnames(reconstruction_world_M$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_N = readRDS('1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_N.rds')
reconstruction_world_N$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_N$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_N$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_N$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_N$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_N$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_N$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_N$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_N$possible_snps = colnames(reconstruction_world_N$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF10 = readRDS('1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF10.rds')
reconstruction_world_ORF10$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF10$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF10$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF10$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF10$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF10$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF10$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF10$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF10$possible_snps = colnames(reconstruction_world_ORF10$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF14 = readRDS('1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF14.rds')
reconstruction_world_ORF14$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF14$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF14$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF14$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF14$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF14$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF14$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF14$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF14$possible_snps = colnames(reconstruction_world_ORF14$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF1a = readRDS('1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF1a.rds')
reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF1a$possible_snps = colnames(reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF1b = readRDS('1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF1b.rds')
reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF1b$possible_snps = colnames(reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF3a = readRDS('1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF3a.rds')
reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF3a$possible_snps = colnames(reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF6 = readRDS('1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF6.rds')
reconstruction_world_ORF6$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF6$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF6$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF6$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF6$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF6$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF6$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF6$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF6$possible_snps = colnames(reconstruction_world_ORF6$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF7a = readRDS('1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF7a.rds')
reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF7a$possible_snps = colnames(reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF7b = readRDS('1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF7b.rds')
reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF7b$possible_snps = colnames(reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF8 = readRDS('1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF8.rds')
reconstruction_world_ORF8$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF8$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF8$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF8$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF8$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF8$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF8$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF8$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF8$possible_snps = colnames(reconstruction_world_ORF8$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF9b = readRDS('1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF9b.rds')
reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF9b$possible_snps = colnames(reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_S = readRDS('1_SARS-CoV-2/2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_S.rds')
reconstruction_world_S$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_S$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_S$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_S$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_S$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_S$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_S$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_S$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_S$possible_snps = colnames(reconstruction_world_S$dataset_with_inferred_reconstruction_codon)[-(1:4)]
########################################################################################################################################

########################################################################################################################################
## For each SNP, look at defining mutation of each group
########################################################################################################################################
association_scores_per_group = function(dataset_with_nodes, dataset_with_inferred_reconstruction, tree, 
                                        possible_snps, upstream_window, downstream_window){
  ## Set list to store results
  group_names = levels(as.factor(dataset_with_nodes$groups))
  scores = as.list(rep(NA, length(group_names)-1))
  n_tips = length(tree$tip.label)
  
  ## For each group (except the initial group, which is the root), look at snp association
  for(j in 1:(length(group_names)-1)){ 
    print(j)
    
    ## Find members of the group and MRCA
    members = dataset_with_nodes$ID[which(dataset_with_nodes$groups == group_names[j])]
    mrca = getMRCA(tree, dataset_with_nodes$name_seq[which(dataset_with_nodes$groups == group_names[j] & dataset_with_nodes$is.node == 'no')])
    members = unique(c(members, mrca))
    
    ## Update members list, with strains downstream, within the time window
    time_mrca = dataset_with_nodes$time[which(dataset_with_nodes$ID == mrca)] ## Reference time for window
    # tmp = getDescendants(tree, mrca)
    # tmp = tmp[which(dataset_with_nodes$time[tmp] < time_mrca + downstream_window)]
    # members = unique(c(tmp, members))
    
    ## Update members list, with strains upstream, within the time window
    ancest = Ancestors(x = tree, node = mrca, type = 'all')
    time_ancest = dataset_with_nodes$time[ancest]
    time_ancest[1] = time_mrca ## Always keep the first ancestor node
    tmp = which(time_ancest < time_mrca - upstream_window)
    if(length(tmp) > 0) ancest = ancest[-tmp]  ## Remove nodes that are outside the time window
    if(length(ancest) == 0) ancest = Ancestors(x = tree, node = mrca, type = 'all')[1]
    groups_ancest = dataset_with_nodes$groups[ancest]
    gr = min(as.numeric(as.character(groups_ancest)))
    ancest = ancest[which(groups_ancest == gr)] ## Makes sure we only have the directly ancestral group, not more
    time_ancest = dataset_with_nodes$time[ancest]
    
    ## Branches of interest
    branches = tree$edge
    branches = branches[match(c(ancest, members), branches[,2]),] ## Take all tips, from ancests, members
    tmp = which(is.na(match(branches[,1], c(ancest, members))))
    if(length(tmp) > 0){  branches = branches[-tmp, ]} ## Remove nodes that are not in ancest or members
    
    ## Branches time
    branches_time = branches
    branches_time[,1] = dataset_with_nodes$time[branches[,1]]
    branches_time[,2] = dataset_with_nodes$time[branches[,2]]
    
    ## Branches group (checked: ok)
    branches_group = branches
    branches_group[,1] = dataset_with_nodes$groups[branches[,1]]
    branches_group[,2] = dataset_with_nodes$groups[branches[,2]]
    
    ## Make branch group matrix binary: 1=group of interest, 0=other group (eg ancestral)
    branches_group[which(branches_group == j, arr.ind = T)] = 1
    branches_group[which(branches_group > j, arr.ind = T)] = 0
    
    snps = time_diff = snps_props_within = snps_props_whole = names_possibles_snps = NULL
    
    for(i in 1:length(possible_snps)){
      branches_snp = branches
      branches_snp[,1] = dataset_with_inferred_reconstruction[branches[,1],
                                                              which(colnames(dataset_with_inferred_reconstruction) == possible_snps[i])]
      branches_snp[,2] = dataset_with_inferred_reconstruction[branches[,2],
                                                              which(colnames(dataset_with_inferred_reconstruction) == possible_snps[i])]
      
      ancestral_state = branches_snp[which.min(branches[,1]),1]
      
      k=1
      ancestral_state = branches_snp[which(branches[,1] == rev(ancest)[k]),1]
      k=2
      while(str_detect(ancestral_state, pattern = 'n|X') == T & k <= length(ancest)-1){
        ancestral_state = branches_snp[which(branches[,1] == rev(ancest)[k]),1]
        k=k+1
      }
      if(str_detect(ancestral_state, pattern = 'n|X') == T){
        ancestral_state = branches_snp[which(branches[,1] == mrca),1][1]
      }
      
      branches_snp[str_detect(branches_snp, pattern = 'n|X')] = ancestral_state
      
      branches_snp_bin = (branches_snp!=ancestral_state)
      
      tmp = which(branches_group[,1] == 1 & branches_group[,2] == 1)
      Px = as.numeric(branches_group[tmp,])
      Sx = as.numeric(branches_snp_bin[tmp,])
      score3 <- (sum((1 - Px)*(1 - Sx), na.rm=TRUE) + sum(Px*Sx, na.rm=TRUE)) / length(Px)
      
      scores[[j]][i] = score3
      
      t = table(branches_snp)
      t = t[which(names(t) != ancestral_state)]
      t = t[which.max(t)]
      names_possibles_snps = c(names_possibles_snps, paste0(ancestral_state, '|', possible_snps[i], '|', names(t)))
    }
    scores[[j]] = scores[[j]]#- median(scores[[j]])
    names(scores[[j]]) = names_possibles_snps
  }
  return(scores)
}
## Time window
upstream_window = 2 # years (not going to the root)
downstream_window = 10 # years (considering all the group)

scores_world_E_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                      dataset_with_inferred_reconstruction = reconstruction_world_E$dataset_with_inferred_reconstruction_codon, 
                                                      tree = tree_sars_world, 
                                                      possible_snps = reconstruction_world_E$possible_snps, 
                                                      upstream_window = upstream_window, 
                                                      downstream_window = downstream_window)
scores_world_M_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                      dataset_with_inferred_reconstruction = reconstruction_world_M$dataset_with_inferred_reconstruction_codon, 
                                                      tree = tree_sars_world, 
                                                      possible_snps = reconstruction_world_M$possible_snps, 
                                                      upstream_window = upstream_window, 
                                                      downstream_window = downstream_window)
scores_world_N_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                      dataset_with_inferred_reconstruction = reconstruction_world_N$dataset_with_inferred_reconstruction_codon, 
                                                      tree = tree_sars_world, 
                                                      possible_snps = reconstruction_world_N$possible_snps, 
                                                      upstream_window = upstream_window, 
                                                      downstream_window = downstream_window)
scores_world_ORF10_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF10$dataset_with_inferred_reconstruction_codon, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF10$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_ORF14_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF14$dataset_with_inferred_reconstruction_codon, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF14$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_ORF1a_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_codon, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF1a$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_ORF1b_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_codon, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF1b$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_ORF3a_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_codon, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF3a$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_ORF6_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                         dataset_with_inferred_reconstruction = reconstruction_world_ORF6$dataset_with_inferred_reconstruction_codon, 
                                                         tree = tree_sars_world, 
                                                         possible_snps = reconstruction_world_ORF6$possible_snps, 
                                                         upstream_window = upstream_window, 
                                                         downstream_window = downstream_window)
scores_world_ORF7a_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_codon, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF7a$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_ORF7b_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_codon, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF7b$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_ORF8_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                         dataset_with_inferred_reconstruction = reconstruction_world_ORF8$dataset_with_inferred_reconstruction_codon, 
                                                         tree = tree_sars_world, 
                                                         possible_snps = reconstruction_world_ORF8$possible_snps, 
                                                         upstream_window = upstream_window, 
                                                         downstream_window = downstream_window)
scores_world_ORF9b_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_codon, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF9b$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_S_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                      dataset_with_inferred_reconstruction = reconstruction_world_S$dataset_with_inferred_reconstruction_codon, 
                                                      tree = tree_sars_world, 
                                                      possible_snps = reconstruction_world_S$possible_snps, 
                                                      upstream_window = upstream_window, 
                                                      downstream_window = downstream_window)

scores_world_E_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                      dataset_with_inferred_reconstruction = reconstruction_world_E$dataset_with_inferred_reconstruction_AA, 
                                                      tree = tree_sars_world, 
                                                      possible_snps = reconstruction_world_E$possible_snps, 
                                                      upstream_window = upstream_window, 
                                                      downstream_window = downstream_window)
scores_world_M_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                      dataset_with_inferred_reconstruction = reconstruction_world_M$dataset_with_inferred_reconstruction_AA, 
                                                      tree = tree_sars_world, 
                                                      possible_snps = reconstruction_world_M$possible_snps, 
                                                      upstream_window = upstream_window, 
                                                      downstream_window = downstream_window)
scores_world_N_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                      dataset_with_inferred_reconstruction = reconstruction_world_N$dataset_with_inferred_reconstruction_AA, 
                                                      tree = tree_sars_world, 
                                                      possible_snps = reconstruction_world_N$possible_snps, 
                                                      upstream_window = upstream_window, 
                                                      downstream_window = downstream_window)
scores_world_ORF10_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF10$dataset_with_inferred_reconstruction_AA, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF10$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_ORF14_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF14$dataset_with_inferred_reconstruction_AA, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF14$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_ORF1a_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_AA, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF1a$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_ORF1b_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_AA, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF1b$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_ORF3a_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_AA, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF3a$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_ORF6_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                         dataset_with_inferred_reconstruction = reconstruction_world_ORF6$dataset_with_inferred_reconstruction_AA, 
                                                         tree = tree_sars_world, 
                                                         possible_snps = reconstruction_world_ORF6$possible_snps, 
                                                         upstream_window = upstream_window, 
                                                         downstream_window = downstream_window)
scores_world_ORF7a_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_AA, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF7a$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_ORF7b_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_AA, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF7b$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_ORF8_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                         dataset_with_inferred_reconstruction = reconstruction_world_ORF8$dataset_with_inferred_reconstruction_AA, 
                                                         tree = tree_sars_world, 
                                                         possible_snps = reconstruction_world_ORF8$possible_snps, 
                                                         upstream_window = upstream_window, 
                                                         downstream_window = downstream_window)
scores_world_ORF9b_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                          dataset_with_inferred_reconstruction = reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_AA, 
                                                          tree = tree_sars_world, 
                                                          possible_snps = reconstruction_world_ORF9b$possible_snps, 
                                                          upstream_window = upstream_window, 
                                                          downstream_window = downstream_window)
scores_world_S_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                      dataset_with_inferred_reconstruction = reconstruction_world_S$dataset_with_inferred_reconstruction_AA, 
                                                      tree = tree_sars_world, 
                                                      possible_snps = reconstruction_world_S$possible_snps, 
                                                      upstream_window = upstream_window, 
                                                      downstream_window = downstream_window)
########################################################################################################################################

save.image('1_SARS-CoV-2/2_analysis_index/3_snp_association/Raw_scores_detection.Rdata')

########################################################################################################################################
## Find significant snps
########################################################################################################################################
edge_lineage_tree = split_world$lineage_tree$edge
edge_lineage_tree[,1] = split_world$tip_and_nodes_groups[match(edge_lineage_tree[,1],names(split_world$tip_and_nodes_groups))]
edge_lineage_tree[,2] = split_world$tip_and_nodes_groups[match(edge_lineage_tree[,2],names(split_world$tip_and_nodes_groups))]
edge_lineage_tree_snps = as.list(edge_lineage_tree)

## Combine all codons
scores_world_sig_all_codons = scores_world_S_codons
sig_threshold = 0.8
for(i in 1:length(scores_world_S_codons)){
  tmp = as.numeric(edge_lineage_tree[which(edge_lineage_tree[,2] == i),])
  if(tmp[1] == length(scores_world_S_codons)+1) {
    ans = NULL
  }else{
    ans = c(paste0('E:', names(which(scores_world_E_codons[[tmp[1]]] > sig_threshold))),
            paste0('M:', names(which(scores_world_M_codons[[tmp[1]]] > sig_threshold))),
            paste0('N:', names(which(scores_world_N_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF10:', names(which(scores_world_ORF10_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF14:', names(which(scores_world_ORF14_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF1a:', names(which(scores_world_ORF1a_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF1b:', names(which(scores_world_ORF1b_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF3a:', names(which(scores_world_ORF3a_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF6:', names(which(scores_world_ORF6_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF7a:', names(which(scores_world_ORF7a_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF7b:', names(which(scores_world_ORF7b_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF8:', names(which(scores_world_ORF8_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF9b:', names(which(scores_world_ORF9b_codons[[tmp[1]]] > sig_threshold))),
            paste0('S:', names(which(scores_world_S_codons[[tmp[1]]] > sig_threshold))))
  }
  des = c(paste0('E:', names(which(scores_world_E_codons[[tmp[2]]] > sig_threshold))),
          paste0('M:', names(which(scores_world_M_codons[[tmp[2]]] > sig_threshold))),
          paste0('N:', names(which(scores_world_N_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF10:', names(which(scores_world_ORF10_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF14:', names(which(scores_world_ORF14_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF1a:', names(which(scores_world_ORF1a_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF1b:', names(which(scores_world_ORF1b_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF3a:', names(which(scores_world_ORF3a_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF6:', names(which(scores_world_ORF6_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF7a:', names(which(scores_world_ORF7a_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF7b:', names(which(scores_world_ORF7b_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF8:', names(which(scores_world_ORF8_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF9b:', names(which(scores_world_ORF9b_codons[[tmp[2]]] > sig_threshold))),
          paste0('S:', names(which(scores_world_S_codons[[tmp[2]]] > sig_threshold))))
  res = des[which(is.na(match(des, ans)))]
  res = res[which(is.na(match(res, c("E:", "M:", "N:", 'ORF10:', 'ORF14:',
                                     'ORF1a:', 'ORF1b:', 'ORF3a:', 'ORF6:', 
                                     'ORF7a:', 'ORF7b:', 'ORF8:', 
                                     'ORF9b:','S:'))))]
  scores_world_sig_all_codons[[i]] = res
}

## Combine all AA
scores_world_sig_all_AA = scores_world_S_AA
sig_threshold = 0.8
for(i in 1:length(scores_world_S_AA)){
  tmp = as.numeric(edge_lineage_tree[which(edge_lineage_tree[,2] == i),])
  if(tmp[1] == length(scores_world_S_AA)+1) {
    ans = NULL
  }else{
    ans = c(paste0('E:', names(which(scores_world_E_AA[[tmp[1]]] > sig_threshold))),
            paste0('M:', names(which(scores_world_M_AA[[tmp[1]]] > sig_threshold))),
            paste0('N:', names(which(scores_world_N_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF10:', names(which(scores_world_ORF10_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF14:', names(which(scores_world_ORF14_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF1a:', names(which(scores_world_ORF1a_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF1b:', names(which(scores_world_ORF1b_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF3a:', names(which(scores_world_ORF3a_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF6:', names(which(scores_world_ORF6_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF7a:', names(which(scores_world_ORF7a_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF7b:', names(which(scores_world_ORF7b_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF8:', names(which(scores_world_ORF8_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF9b:', names(which(scores_world_ORF9b_AA[[tmp[1]]] > sig_threshold))),
            paste0('S:', names(which(scores_world_S_AA[[tmp[1]]] > sig_threshold))))
  }
  des = c(paste0('E:', names(which(scores_world_E_AA[[tmp[2]]] > sig_threshold))),
          paste0('M:', names(which(scores_world_M_AA[[tmp[2]]] > sig_threshold))),
          paste0('N:', names(which(scores_world_N_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF10:', names(which(scores_world_ORF10_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF14:', names(which(scores_world_ORF14_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF1a:', names(which(scores_world_ORF1a_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF1b:', names(which(scores_world_ORF1b_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF3a:', names(which(scores_world_ORF3a_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF6:', names(which(scores_world_ORF6_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF7a:', names(which(scores_world_ORF7a_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF7b:', names(which(scores_world_ORF7b_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF8:', names(which(scores_world_ORF8_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF9b:', names(which(scores_world_ORF9b_AA[[tmp[2]]] > sig_threshold))),
          paste0('S:', names(which(scores_world_S_AA[[tmp[2]]] > sig_threshold))))
  res = des[which(is.na(match(des, ans)))]
  res = res[which(is.na(match(res, c("E:", "M:", "N:", 'ORF10:', 'ORF14:',
                                     'ORF1a:', 'ORF1b:', 'ORF3a:', 'ORF6:', 
                                     'ORF7a:', 'ORF7b:', 'ORF8:', 
                                     'ORF9b:','S:'))))]
  scores_world_sig_all_AA[[i]] = res
}

## Get match with Obermeyer
data_obermeyer = read.csv('1_SARS-CoV-2/1_2022_Obermeyer_mutations/mutations.tsv', sep = '\t')
data_obermeyer$gene = unlist(lapply(data_obermeyer$mutation, function(x)str_split(x, ':')[[1]][1]))
data_obermeyer$mutation = str_replace_all(data_obermeyer$mutation, pattern = 'STOP', replacement = '*')
tmp = unique(toupper(unlist(scores_world_sig_all_AA)))
tmp = str_remove_all(tmp, '\\|')           
m = match(tmp, toupper(data_obermeyer$mutation))
length(which(is.na(m)))/length(m)
1-length(which(is.na(m)))/length(m)
length(m)
data_obermeyer$found = rep(1, nrow(data_obermeyer))
data_obermeyer$found[m] = 2
data_obermeyer_match = data_obermeyer[m,]
########################################################################################################################################

########################################################################################################################################
## Plot comparison Obermeyer's list
########################################################################################################################################
plot_Obermeyer_match = function(){
  par(oma = c(0,0,0,0), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
      cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)
  plot(NULL, xlim = c(0,3000), ylim = c(-0.2, 0.35), 
       bty = 'n', 
       yaxt = 'n',
       main = 'All mutations',
       xlab = 'Ranked mutations',
       ylab = 'Δ log R')
  axis(2, las = 2, mgp = c(1.5,0.5,0))
  idx = which(data_obermeyer$found == 1)
  points(idx, data_obermeyer$Δ.log.R[idx], col = 'black', pch = 16, cex = 0.5)
  idx = which(data_obermeyer$found == 2)
  points(idx, data_obermeyer$Δ.log.R[idx], col = 'firebrick', pch = 16, cex = 0.5)
  
  legend('topright', pch = 16, col = c('firebrick', 'black'), legend = c('yes', 'no'), title = 'Found in automatic detection', bty = 'n', 
         cex = 0.7)
}
plot_Obermeyer_match_zoom = function(){
  par(oma = c(0,0,0,0), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
      cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)
  plot(NULL, xlim = c(0,200), ylim = c(0, 0.4), 
       bty = 'n', 
       main = 'Top ranked mutations',yaxt = 'n',
       xlab = 'Ranked mutations',
       ylab = 'Δ log R')
  axis(2, las = 2, mgp = c(1.5,0.5,0))
  idx = which(data_obermeyer$found == 1)
  points(idx, data_obermeyer$Δ.log.R[idx], col = 'black', pch = 16, cex = 0.35)
  idx = which(data_obermeyer$found == 2)
  points(idx, data_obermeyer$Δ.log.R[idx], col = 'firebrick', pch = 16, cex = 0.35)
  axis(2, las = 2, mgp = c(1.5,0.5,0))
}
ggdraw(plot_Obermeyer_match_zoom)
ggdraw(plot_Obermeyer_match)
########################################################################################################################################

########################################################################################################################################
## Full genome, substitutions (NS, ie AA change) with signal
########################################################################################################################################
genes = read.csv('1_SARS-CoV-2/1_refseq/Position_ORFs_nextstrain.csv')
genes = rbind(genes, 
              c(15, NA, NA, 'region', 21563+(333-1)*3, 21563+(527-1)*3+2, NA, '+', NA, 'region_name "RBD"', "RBD"))
genes_to_plot = genes
genes_to_plot = genes_to_plot[-which(genes_to_plot$gene == 'ORF3a' | genes_to_plot$gene == 'ORF10' | genes_to_plot$gene == 'ORF14'| genes_to_plot$gene == 'ORF6' | genes_to_plot$gene == 'ORF7a'|
                                       genes_to_plot$gene == 'ORF7b' | genes_to_plot$gene == 'ORF8'| genes_to_plot$gene == 'ORF9b'),]
colors = MetBrewer::met.brewer(name="Derain", n=15)

scores_world_full_AA = NULL
for(i in 1:(nrow(genes)-1)){
  if(genes$gene[i] == 'E') {
    tmp = apply(do.call(rbind, scores_world_E_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'M') {
    tmp = apply(do.call(rbind, scores_world_M_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'N') {
    tmp = apply(do.call(rbind, scores_world_N_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF10') {
    tmp = apply(do.call(rbind, scores_world_ORF10_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF14') {
    tmp = apply(do.call(rbind, scores_world_ORF14_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF1a') {
    tmp = apply(do.call(rbind, scores_world_ORF1a_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF1b') {
    tmp = apply(do.call(rbind, scores_world_ORF1b_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF3a') {
    tmp = apply(do.call(rbind, scores_world_ORF3a_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF6') {
    tmp = apply(do.call(rbind, scores_world_ORF6_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF7a') {
    tmp = apply(do.call(rbind, scores_world_ORF7a_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF7b') {
    tmp = apply(do.call(rbind, scores_world_ORF7b_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF8') {
    tmp = apply(do.call(rbind, scores_world_ORF8_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF9b') {
    tmp = apply(do.call(rbind, scores_world_ORF9b_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'S') {
    tmp = apply(do.call(rbind, scores_world_S_AA), MAR = 2, max)
  } 
  names_tmp = as.numeric(unlist(lapply(names(tmp), function(x)str_split(x, '\\|')[[1]][2])))
  names_tmp = names_tmp*3-1.5 + as.numeric(genes$start[i]) - 1
  names(tmp) = names_tmp
  scores_world_full_AA = c(scores_world_full_AA, tmp)
}

plot_scores_genome_density_AA = function(){
  par(oma = c(0,0,0,0), mar = c(2,2,0,0), mgp = c(0,0.1,0), family = 'Arial',
      cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)
  dens = density(as.numeric(names(scores_world_full_AA)), weights = scores_world_full_AA/sum(scores_world_full_AA), 
                 bw = 50, n = 3000)
  plot(NULL, xlim = c(0, 30000), ylim = c(0, max(dens$y)), bty = 'n', yaxt = 'n',
       ylab = '')
  axis(2, las = 2)
  for(i in 1:nrow(genes)){
    if(!is.na(match(genes$gene[i], genes_to_plot$gene))){
      polygon(x = c(c(genes$start[i], genes$end[i]), rev(c(genes$start[i], genes$end[i]))), y = c(0, 0, 1, 1), 
              border = F, col = adjustcolor(colors[i], alpha.f = 0.25))
    }
  }
  polygon(x = c(dens$x, rev(dens$x)), 
          y = c(dens$y, rep(0, length(dens$y))), 
          border = F, col = 'grey30')
}
ggdraw(plot_scores_genome_density_AA)
############################################

########################################################################################################################################
## Create meaningful lineage trees
########################################################################################################################################
create_timed_lineage_tree = function(condensed_tree, node_group, dataset, time_root = 1900, tip_date = NULL){
  ## TO DO: clean code
  if(length(which(is.na(node_group))) > 0){
    tip_to_remove = which(is.na(condensed_tree$tip.label))
    condensed_tree$tip.label = condensed_tree$tip.label[-tip_to_remove]
    condensed_tree$edge = condensed_tree$edge[-which(condensed_tree$edge[,2] == tip_to_remove),]
    condensed_tree$edge[which(condensed_tree$edge<= length(condensed_tree$tip.label), arr.ind = T)] = condensed_tree$edge[which(condensed_tree$edge<= length(condensed_tree$tip.label), arr.ind = T)]-1
    condensed_tree$edge[which(condensed_tree$edge> length(condensed_tree$tip.label), arr.ind = T)] = condensed_tree$edge[which(condensed_tree$edge> length(condensed_tree$tip.label), arr.ind = T)]-1
    # condensed_tree = drop.tip(condensed_tree, tip = which(is.na(node_group)))
    node_group = node_group[-tip_to_remove]
  }
  
  tree_to_time = condensed_tree
  edge_group = condensed_tree$edge
  edge_group[,1] = as.numeric(node_group[condensed_tree$edge[,1]])
  edge_group[,2] = as.numeric(node_group[condensed_tree$edge[,2]])
  
  ####################################################
  ## Step 1: find MRCA each group and reorder edges to go fro oldest to newest
  ####################################################
  t_mrca_group = NULL
  groups_names = unique(node_group)
  for(j in 1:length(groups_names)){
    t_mrca_group = c(t_mrca_group, min(dataset$time[which(dataset$groups == groups_names[j])]))
  }
  names(t_mrca_group) = groups_names
  
  o = order(t_mrca_group)
  tmp = as.numeric(t_mrca_group[match(edge_group[,2], names(t_mrca_group))])
  tree_to_time$edge = tree_to_time$edge[order(tmp),]
  edge_group = edge_group[order(tmp),]
  ####################################################
  
  ####################################################
  ## Step 2: add tips for groups that are not tips
  ####################################################
  Ntips = length(tree_to_time$tip.label) ## Warning: tree not binary (yet)!
  Nnodes = tree_to_time$Nnode
  New_N_tips = Ntips + tree_to_time$Nnode
  tree_to_time$edge[which(tree_to_time$edge > Ntips, arr.ind = T)] = tree_to_time$edge[which(tree_to_time$edge > Ntips, arr.ind = T)] - Ntips  + New_N_tips
  edge_group_new = edge_group
  k = max(tree_to_time$edge) + 1
  ## Start adding a node at the root
  tree_to_time$edge = rbind(c(tree_to_time$edge[1,1], Ntips + 1), tree_to_time$edge)
  edge_group_new = rbind(c(edge_group_new[1,1], edge_group_new[1,1]), edge_group_new)
  
  for(i in 2:Nnodes){
    ## Looking for nodes in the tree
    index = which(tree_to_time$edge[,2] == New_N_tips + i)
    number_nodes = length(which(tree_to_time$edge[,1] == New_N_tips + i))
    print(number_nodes)
    ## Adding a tip for this node
    if(index+1 <= nrow(tree_to_time$edge)){
      tree_to_time$edge = rbind(tree_to_time$edge[1:index,],
                                c(tree_to_time$edge[index,2], Ntips + i),
                                tree_to_time$edge[(index+1):nrow(tree_to_time$edge),])
      edge_group_new = rbind(edge_group_new[1:index,], 
                             c(edge_group_new[index,2], edge_group_new[index,2]),
                             edge_group_new[(index+1):nrow(edge_group_new),])
    }else{
      tree_to_time$edge = rbind(tree_to_time$edge[1:index,],
                                c(tree_to_time$edge[index,2], Ntips + i))
      edge_group_new = rbind(edge_group_new[1:index,], 
                             c(edge_group_new[index,2], edge_group_new[index,2]))
    }
  }
  ## Make it a phylo object
  tips = as.numeric(names(table(tree_to_time$edge))[which(table(tree_to_time$edge) == 1)])
  tree_to_time$tip.label = edge_group_new[match(tips, tree_to_time$edge[,2]),2]
  tree_to_time$edge[match(tips, tree_to_time$edge)] =  1:length(tips)
  nodes = unique(tree_to_time$edge[which(tree_to_time$edge > length(tips))])
  tree_to_time$Nnode = length(nodes)
  edge_tmp = tree_to_time$edge
  for(i in 1:tree_to_time$Nnode){ ## Renumber nodes, in ascending order
    edge_tmp[which(!is.na(match(tree_to_time$edge, nodes[i])))] = as.integer(length(tips)+i)
  }
  tree_to_time$edge = edge_tmp
  attributes(tree_to_time) = list('class' = 'phylo', 'names' = names(tree_to_time))
  tree_to_time = as.phylo(tree_to_time) ## Create phylo object
  
  ## Output group of each  node
  node_groups = as.character(edge_group_new[match(1:(length(tree_to_time$tip.label)+tree_to_time$Nnode), tree_to_time$edge)])
  ####################################################
  
  ####################################################
  ## Step 3: Make tree non-binary and timed
  ####################################################
  timed_tree = tree_to_time
  edge_group_timed_tree = edge_group_new
  ## First, make the tree non-binary
  ## If root is non-binary: 
  if(length(which(timed_tree$edge[,1] == length(timed_tree$tip.label)+1)) > 2){
    id = which(timed_tree$edge[,1] == length(timed_tree$tip.label)+1)
    l = length(which(timed_tree$edge == length(timed_tree$tip.label)+1))
    k = 0.001
    for(j in 3:l){
      ## For those groups, update name of splitting node
      timed_tree$edge[id[j],1] = timed_tree$edge[id[1],1]+k
      k = k + 0.001
    }
    
    ## Then split the main branch into as many pieces as needed
    branches_to_add = rbind(c(timed_tree$edge[id[1],1],timed_tree$edge[id[1],1]+0.001))
    groups_to_add = rbind(c(edge_group_timed_tree[id[1],1],edge_group_timed_tree[id[1],1]))
    k = 0.001
    if(l > 3){
      for(j in 4:l){
        branches_to_add = rbind(branches_to_add,
                                c(timed_tree$edge[id[1],1]+k, timed_tree$edge[id[1],1]+k+0.001))
        groups_to_add = rbind(groups_to_add,
                              c(edge_group_timed_tree[id[1],1],edge_group_timed_tree[id[1],1]))
        k = k + 0.001
      }
    }
    branches_to_add = rbind(branches_to_add,
                            c(timed_tree$edge[id[1],1]+k, timed_tree$edge[id[1],2]))
    groups_to_add = rbind(groups_to_add,
                          c(edge_group_timed_tree[id[1],1],edge_group_timed_tree[id[1],2]))
    
    if(id[1] != 1){
      timed_tree$edge = rbind(timed_tree$edge[1:(id[1]-1),],
                              branches_to_add,
                              timed_tree$edge[(id[1]+1):nrow(timed_tree$edge),])
      edge_group_timed_tree = rbind(edge_group_timed_tree[1:(id[1]-1),],
                                    groups_to_add,
                                    edge_group_timed_tree[(id[1]+1):nrow(edge_group_timed_tree),])
    }else{
      timed_tree$edge = rbind(branches_to_add,
                              timed_tree$edge[(id[1]+1):nrow(timed_tree$edge),])
      edge_group_timed_tree = rbind(groups_to_add,
                                    edge_group_timed_tree[(id[1]+1):nrow(edge_group_timed_tree),])
    }
    
  }
  ## Other nodes
  tmp = names(which(table(timed_tree$edge) > 3))
  if(length(tmp) > 0){
    for(i in 1:length(tmp)){
      id = which(timed_tree$edge[,1] == tmp[i])
      time_id = t_mrca_group[match(edge_group_timed_tree[id,2],names(t_mrca_group))]
      id = id[order(time_id)]
      time_id = time_id[order(time_id)]
      l = length(id)
      
      k = 0.001
      for(j in 3:l){
        ## For those groups, update name of splitting node
        timed_tree$edge[id[j],1] = timed_tree$edge[id[1],1]+k
        k = k + 0.001
      }
      
      ## Then split the main branch into as many pieces as needed
      branches_to_add = rbind(c(timed_tree$edge[id[1],1],timed_tree$edge[id[1],1]+0.001))
      groups_to_add = rbind(c(edge_group_timed_tree[id[1],1],edge_group_timed_tree[id[1],1]))
      k = 0.001
      if(l > 3){
        for(j in 4:l){
          branches_to_add = rbind(branches_to_add, 
                                  c(timed_tree$edge[id[1],1]+k, timed_tree$edge[id[1],1]+k+0.001))
          groups_to_add = rbind(groups_to_add,
                                c(edge_group_timed_tree[id[1],1],edge_group_timed_tree[id[1],1]))
          k = k + 0.001
        }
      }
      branches_to_add = rbind(branches_to_add, 
                              c(timed_tree$edge[id[1],1]+k, timed_tree$edge[id[1],2]))
      groups_to_add = rbind(groups_to_add,
                            c(edge_group_timed_tree[id[1],1],edge_group_timed_tree[id[1],2]))
      
      timed_tree$edge = rbind(timed_tree$edge[1:(id[1]-1),],
                              branches_to_add,
                              timed_tree$edge[(id[1]+1):nrow(timed_tree$edge),])
      edge_group_timed_tree = rbind(edge_group_timed_tree[1:(id[1]-1),],
                                    groups_to_add,
                                    edge_group_timed_tree[(id[1]+1):nrow(edge_group_timed_tree),])
    }
  }
  
  ## Make it a phylo object
  tips = as.numeric(names(table(timed_tree$edge))[which(table(timed_tree$edge) == 1)])
  timed_tree$tip.label = edge_group_timed_tree[match(tips, timed_tree$edge[,2]),2]
  timed_tree$edge[match(tips, timed_tree$edge)] =  1:length(tips)
  nodes = unique(as.character(timed_tree$edge[which(timed_tree$edge > length(tips))]))
  timed_tree$Nnode = length(nodes)
  edge_tmp = timed_tree$edge
  for(i in 1:timed_tree$Nnode){ ## Renumber nodes, in ascending order
    edge_tmp[which(!is.na(match(timed_tree$edge, nodes[i])))] = as.integer(length(tips)+i)
  }
  timed_tree$edge = edge_tmp
  attributes(timed_tree) = list('class' = 'phylo', 'names' = names(timed_tree))
  timed_tree = as.phylo(timed_tree) ## Create phylo object
  
  ## Second, compute edge lengths
  ## Compute time of each node
  Ntips = timed_tree$Nnode + 1
  list_times = NULL
  ## For the root
  # list_times = c(list_times, time_root)
  for(i in 1:timed_tree$Nnode){
    index = which(timed_tree$edge[,1] == i+Ntips)
    j = which(edge_group_timed_tree[index,1] - edge_group_timed_tree[index,2] != 0)
    if(length(j) > 1) print('problem')
    if(edge_group_timed_tree[index[j],1] != edge_group_timed_tree[index[j],2]){
      list_times = c(list_times, t_mrca_group[which(groups_names == edge_group_timed_tree[index[j],2])])
    }
  }
  ## Replace these times in the edge matrix
  edge_dates = matrix(NA, nrow = nrow(timed_tree$edge), ncol = 2)
  edge_dates[,1] = list_times[match(timed_tree$edge[,1], (1:timed_tree$Nnode)+Ntips)]
  edge_dates[,2] = list_times[match(timed_tree$edge[,2], (1:timed_tree$Nnode)+Ntips)]
  ## Add times of the tips
  index_tips = which(timed_tree$edge <= length(timed_tree$tip.label), arr.ind = T)
  if(is.null(tip_date)){
    for(ind in index_tips[,1]){
      gr = edge_group_timed_tree[ind,2]
      edge_dates[ind,2] = max(dataset$time[which(dataset$groups == gr)])
    }
  }else{
    edge_dates[index_tips[,1],2] = tip_date # Add dates of tips
  }
  ## Compute edge length and include it in the tree
  edge_length = edge_dates[,2] - edge_dates[,1]
  timed_tree$edge.length = edge_length
  node_groups = as.character(edge_group_timed_tree[match(1:(length(timed_tree$tip.label)+timed_tree$Nnode), timed_tree$edge)])
  names(node_groups) = 1:(length(timed_tree$tip.label)+timed_tree$Nnode)
  
  return(list('lineage_tree' = timed_tree,
              'node_group' = node_groups,
              'edge_group' = edge_group_timed_tree,
              'mrca' = t_mrca_group))
}

lineage_tree_sarscov2_world = create_timed_lineage_tree(split_world$lineage_tree, 
                                                        split_world$tip_and_nodes_groups, 
                                                        dataset_with_nodes_world)
########################################################################################################################################

########################################################################################################################################
## Plot lineage, with AA annotations
########################################################################################################################################
lineage_tree = lineage_tree_sarscov2_world
groups = levels(factor(split_world$tip_and_nodes_groups))
sig_snp_nodes = rep(NA, length(lineage_tree$node_group))
for(i in 1:(length(groups)-1)){
  idx = which(as.character(lineage_tree$node_group) == groups[i])
  tmp = scores_world_sig_all_AA[[i]]
  # tmp = tmp[which(tmp$syn == F | is.na(tmp$syn)),]
  # tmp$gene_name[which(is.na(tmp$gene_name))] = tmp$SNP[which(is.na(tmp$gene_name))]
  if(length(idx) == 1){
    sig_snp_nodes[idx[1]] = paste0(tmp, collapse = '/')
  }else if(length(idx) == 2){
    sig_snp_nodes[idx[2]] = paste0(tmp, collapse = '/')
    sig_snp_nodes[idx[1]] = ''
  }else if(length(idx) > 2){
    sig_snp_nodes[idx[2]] = paste0(tmp, collapse = '/')
    sig_snp_nodes[idx[1]] = ''
    sig_snp_nodes[idx[3:length(idx)]] = ''
  }
}

sig_snp_nodes[which(is.na(sig_snp_nodes))] = ''
names(colors_groups_world) = seq_along(colors_groups_world)
iGroup = seq_along(colors_groups_world)
names(iGroup) = seq_along(colors_groups_world)
p = ggtree(lineage_tree$lineage_tree, aes(col = lineage_tree$node_group), layout="roundrect", size=0.5, mrsd="2023-04-20") +
  geom_nodepoint(size=1.5, alpha=1, aes(col = lineage_tree$node_group)) +
  geom_tippoint(size=1.5, alpha=1, aes(col = lineage_tree$node_group)) +
  geom_text(size=1.5, aes(label = sig_snp_nodes), hjust = 1, vjust = -2)+
  scale_x_continuous(limits = c(2019.75, 2023.5), breaks = seq(2020, 2023.5, 0.5), labels = seq(2020, 2023.5, 0.5))+
  scale_color_manual(values=colors_groups_world, breaks = iGroup, na.value = 'white', name = 'Groups') +
  theme_tree2()
p
########################################################################################################################################

########################################################################################################################################
## Save panel
########################################################################################################################################
panel = plot_grid(p, ggdraw(plot_scores_genome_density_AA),ggdraw(plot_Obermeyer_match),
                  rel_widths = c(1, 1.5, 1.25), labels = '', ncol = 3)

ggsave(filename = '1_SARS-CoV-2/2_analysis_index/3_snp_association/Panel_association_sarscov2_world.pdf', plot = panel, device = 'pdf', scale = 1,
       width = 22, height = 5, units = 'cm')
########################################################################################################################################






########################################################################################################################################
## Plot Spike mutations next to tree
########################################################################################################################################
scores_world_sig_all_AA
gene_to_plot = 'S'
dataset_with_inferred_reconstruction = reconstruction_world_S$dataset_with_inferred_reconstruction_AA

AA_mat_total = matrix(rep(100, length(which(dataset_with_inferred_reconstruction$is.node == 'no'))), ncol = 1)

for(group_to_plot in 1:13){
  # group_to_plot = 1
  
  AAs = scores_world_sig_all_AA[[group_to_plot]]
  AAs = AAs[which(unlist(lapply(AAs, function(x)str_split(x, ':')[[1]][1])) == gene_to_plot)]
  if(length(AAs) >= 1){
    AAs_pos = as.numeric(unlist(lapply(AAs, function(x)str_split(x, '\\|')[[1]][2])))
    AAs_ref = unlist(lapply(AAs, function(x)str_split(str_split(x, '\\|')[[1]][1], ':')[[1]][2]))
    AAs_var = unlist(lapply(AAs, function(x)str_split(x, '\\|')[[1]][3]))
    
    AA_mat = matrix(dataset_with_inferred_reconstruction[which(dataset_with_inferred_reconstruction$is.node == 'no'),
                                                         which(colnames(dataset_with_inferred_reconstruction) == AAs_pos[1])], ncol = 1)
    AA_mat[which(AA_mat[,1] != AAs_ref[1] & AA_mat[,1] != AAs_var[1]),1] = NA
    AA_mat[which(AA_mat[,1] == AAs_var[1]),1] = group_to_plot
    AA_mat[which(AA_mat[,1] == AAs_ref[1]),1] = 0
    if(length(AAs) > 1){
      for(i in 2:length(AAs)){
        AA_mat = cbind(AA_mat, dataset_with_inferred_reconstruction[which(dataset_with_inferred_reconstruction$is.node == 'no'),
                                                                    which(colnames(dataset_with_inferred_reconstruction) == AAs_pos[i])])
        AA_mat[which(AA_mat[,i] != AAs_ref[i] & AA_mat[,i] != AAs_var[i]),i] = NA
        AA_mat[which(AA_mat[,i] == AAs_var[i]),i] = group_to_plot
        AA_mat[which(AA_mat[,i] == AAs_ref[i]),i] = 0
      }
    }
    AA_mat = cbind(AA_mat,
                   matrix(rep(100, length(which(dataset_with_inferred_reconstruction$is.node == 'no'))), ncol = 1))
    
    colnames(AA_mat) = c(paste0(group_to_plot, ':', AAs), group_to_plot) 
    rownames(AA_mat) = dataset_with_inferred_reconstruction$name_seq[which(dataset_with_inferred_reconstruction$is.node == 'no')]
    
    AA_mat_total = cbind(AA_mat_total, 
                         AA_mat)
  }
}

cols = colors_groups_world
names(cols) = as.character(1:max(as.numeric(name_groups_world)))

col_AA = c('white', colors_groups_world, 'grey')
names(col_AA) = c(0, as.character(1:max(as.numeric(name_groups_world))), 100)

plot_tree_sars_world_groups <- ggtree(tree_sars_world, mrsd=lubridate::date_decimal(max(times_seqs_world)), size = 0.08,
                                      aes(color = as.character(dataset_with_nodes_world$groups))) + 
  scale_color_manual(values = cols)+
  ggplot2::ylim(-0.5, NA) + 
  coord_cartesian(clip = "off")
  # theme_tree2()

plot_tree_sars_world_groups = gheatmap(plot_tree_sars_world_groups, 
                                       AA_mat_total, 
                                       offset=0.1, width=10, 
                                       colnames=T, 
                                       colnames_position="top",
                                       colnames_angle=90, hjust=0, font.size=2.5,
                                       legend_title="Group", color=NA) +
  scale_fill_manual(values = col_AA, na.value = 'white')+
  coord_cartesian(clip = "off") +
  theme(legend.position = 'none')

ggsave(filename = '1_SARS-CoV-2/2_analysis_index/3_snp_association/Tree_with_S_mutations_sarscov2_world.pdf', plot = plot_tree_sars_world_groups, device = 'pdf', scale = 1,
       width = 27, height = 18, units = 'cm')
########################################################################################################################################


