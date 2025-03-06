## Codes for finding lineage defining mutations, B. pertussis
## Noemie Lefrancq
## Code cleaned 06/03/2025
########################################################################################################################################
## Packages used
########################################################################################################################################
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

########################################################################################################################################
## Useful functions
########################################################################################################################################
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
########################################################################################################################################

########################################################################################################################################
## Load data
########################################################################################################################################
load('3_BordetellaPertussis/2_analysis_index/1_index_computations/Initial_index_computation_and_parameters.Rdata')
load('3_BordetellaPertussis/2_analysis_index/2_find_index_groups/Lineages_detected.Rdata')
########################################################################################################################################


########################################################################################################################################
## PART 1: ancestral state reconstruction: reconstruct the SNP of each. node
########################################################################################################################################
reconstruct_node_states = function(tree, dataset_tips, dataset_with_nodes, min_prop, max_prop, names_seqs){
  # Meta-data with nodes
  dataset_with_inferred_resonstruction = dataset_with_nodes[,1:4]
  possible_snps = names(dataset_tips[,4:ncol(dataset_tips)])

  ## Filter SNP on frequency: only keep SNP that are present in:
  prevalence = t(apply(dataset_tips[,4:ncol(dataset_tips)], MARGIN = 2, function(x)table(factor(x, levels = c("0", "1")))))
  prevalence_prop = prevalence[,2]/(prevalence[,1]+prevalence[,2])
  prevalence = cbind(prevalence, prevalence_prop)
  possible_snps = names(which(prevalence[,3] >= min_prop & prevalence[,3] <= max_prop))

  ## Filter dataset_tips accordingly
  a = which(is.na(match(colnames(dataset_tips), possible_snps)) == F)
  dataset_tips = dataset_tips[,c(1:3, a)]

  ## Check the tree has no node label (otherwise ace produces an error)
  if(!is.null(tree$node.label)){
    tree$node.label = NULL
  }
  for(i in possible_snps){
    print(i)

    snp_data = dataset_tips[,which(colnames(dataset_tips) == i)]
    snp_data = as.factor(snp_data)
    names(snp_data) = names_seqs

    rec = try(ape::ace(phy = tree, x = snp_data, type = "discrete", method = "ML"), silent=TRUE)
    if("try-error" %in% class(rec) | is.null(rec)) {
      rec_all = c(snp_data, rep(NA, tree$Nnode))
      dataset_with_inferred_resonstruction = cbind(dataset_with_inferred_resonstruction,  rec_all)
    }else{
      tmp = rec$lik.anc[,2]
      snp_data = as.numeric(as.character(snp_data))
      rec_all = c(snp_data, tmp) ## List of all states: first all tips, then all nodes

      ## Find first state, to set it to 0, always
      first_state = rec_all[which(dataset_with_inferred_resonstruction$ID == length(names_seqs) + 1)]
      first_state = round(first_state, digits = 0)

      ## Write reconstruction in the big dataset
      if(first_state < 0.5){ ## First state is 0: all good
        dataset_with_inferred_resonstruction = cbind(dataset_with_inferred_resonstruction, rec_all)
      }
      if(first_state > 0.5){ ## First state is 1: have to change to 0
        dataset_with_inferred_resonstruction = cbind(dataset_with_inferred_resonstruction,  1 - rec_all)
      }
    }
    ## Set column name to position name
    colnames(dataset_with_inferred_resonstruction)[which(colnames(dataset_tips) == i) +1] = i
  }
  return(list('dataset_with_inferred_resonstruction' = dataset_with_inferred_resonstruction,
              'snp_prevalence' = prevalence,
              'possible_snps' = possible_snps))
}

reconstruction_France = reconstruct_node_states(tree = tree_France,
                                                dataset_tips = dataset_tips_France,
                                                dataset_with_nodes = dataset_with_nodes_France,
                                                min_prop = 0.01, ## to speed up you can consider less SNPS by putting eg min_prop = 0.05
                                                max_prop = 0.99, ## to speed up you can consider less SNPS by putting eg min_prop = 0.95
                                                names_seqs = names_seqs_France)

## Annotation: For all possible SNP, find the gene name + position of mutation
## Data on Tohama and translation
data_genes$start = (length_genome-data_genes$start+1)
data_genes$end = (length_genome-data_genes$end+1)
translation_table = matrix(c('TTT', 'F', 'Phe', 'TCT', 'S', 'Ser', 'TAT', 'Y', 'Tyr', 'TGT', 'C', 'Cys',
                             'TTC', 'F', 'Phe', 'TCC', 'S', 'Ser', 'TAC', 'Y', 'Tyr', 'TGC', 'C', 'Cys',
                             'TTA', 'L', 'Leu', 'TCA', 'S', 'Ser', 'TAA', '*', 'Ter', 'TGA', '*', 'Ter',
                             'TTG', 'L', 'Leu', 'TCG', 'S', 'Ser', 'TAG', '*', 'Ter', 'TGG', 'W', 'Trp',
                             'CTT', 'L', 'Leu', 'CCT', 'P', 'Pro', 'CAT', 'H', 'His', 'CGT', 'R', 'Arg',
                             'CTC', 'L', 'Leu', 'CCC', 'P', 'Pro', 'CAC', 'H', 'His', 'CGC', 'R', 'Arg',
                             'CTA', 'L', 'Leu', 'CCA', 'P', 'Pro', 'CAA', 'Q', 'Gln', 'CGA', 'R', 'Arg',
                             'CTG', 'L', 'Leu', 'CCG', 'P', 'Pro', 'CAG', 'Q', 'Gln', 'CGG', 'R', 'Arg',
                             'ATT', 'I', 'Ile', 'ACT', 'T', 'Thr', 'AAT', 'N', 'Asn', 'AGT', 'S', 'Ser',
                             'ATC', 'I', 'Ile', 'ACC', 'T', 'Thr', 'AAC', 'N', 'Asn', 'AGC', 'S', 'Ser',
                             'ATA', 'I', 'Ile', 'ACA', 'T', 'Thr', 'AAA', 'K', 'Lys', 'AGA', 'R', 'Arg',
                             'ATG', 'M', 'Met', 'ACG', 'T', 'Thr', 'AAG', 'K', 'Lys', 'AGG', 'R', 'Arg',
                             'GTT', 'V', 'Val', 'GCT', 'A', 'Ala', 'GAT', 'D', 'Asp', 'GGT', 'G', 'Gly',
                             'GTC', 'V', 'Val', 'GCC', 'A', 'Ala', 'GAC', 'D', 'Asp', 'GGC', 'G', 'Gly',
                             'GTA', 'V', 'Val', 'GCA', 'A', 'Ala', 'GAA', 'E', 'Glu', 'GGA', 'G', 'Gly',
                             'GTG', 'V', 'Val', 'GCG', 'A', 'Ala', 'GAG', 'E', 'Glu', 'GGG', 'G', 'Gly'),
                           ncol = 3, byrow = T)

annotate_possible_snps_tohama = function(possible_snps, data_vcf){
  possible_snps_genes = rep(NA, length(possible_snps))
  possible_snps_product = rep(NA, length(possible_snps))
  possible_snps_gene_ID_old = rep(NA, length(possible_snps))
  possible_snps_gene_function = rep(NA, length(possible_snps))
  possible_snps_gene_prod_localisation = rep(NA, length(possible_snps))
  possible_snps_pos_in_gene = rep(NA, length(possible_snps))
  possible_snps_codon_ref = rep(NA, length(possible_snps))
  possible_snps_codon_snp = rep(NA, length(possible_snps))
  possible_snps_close_upstream_gene = rep(NA, length(possible_snps))
  possible_snps_close_upstream_gene_strand = rep(NA, length(possible_snps))
  possible_snps_close_upstream_gene_distance = rep(NA, length(possible_snps))
  possible_snps_close_downstream_gene = rep(NA, length(possible_snps))
  possible_snps_close_downstream_gene_strand = rep(NA, length(possible_snps))
  possible_snps_close_downstream_gene_distance = rep(NA, length(possible_snps))
  for(j in 1:length(possible_snps)){
    ## Code checked for a large random selection of snp, all good (esp. strands + and -)
    pos_name = possible_snps[j]
    pos_numeric = as.numeric(str_split(possible_snps[j], pattern = '_')[[1]][1])
    a = which(data_genes$start>=pos_numeric & data_genes$end<=pos_numeric)
    tmp = NA
    tmp2 = NA
    tmp3 = NA
    tmp4 = NA
    tmp5 = NA
    tmp6 = NA
    if(length(a)>0){ ## TO DO: consider that one position can be in multiple genes (+ and -)
      tmp = data_genes$locus_tag[a]
      tmp2 = data_genes$old_locus_tag[a]
      tmp3 = data_genes$gene[a]
      if(data_genes$strand[a] == '+'){
        tmp4 = data_genes$start[a] - pos_numeric + 1
        if(tmp4%%3 == 1){
          tmp5 = paste0(tohama_full_genome[pos_numeric], tohama_full_genome[pos_numeric - 1], tohama_full_genome[pos_numeric - 2])
          tmp6 = paste0(tolower(data_vcf$ALT[which(data_vcf$POS == pos_name)]), tohama_full_genome[pos_numeric - 1], tohama_full_genome[pos_numeric - 2])
        }
        if(tmp4%%3 == 2){
          tmp5 = paste0(tohama_full_genome[pos_numeric + 1], tohama_full_genome[pos_numeric], tohama_full_genome[pos_numeric - 1])
          tmp6 = paste0(tohama_full_genome[pos_numeric + 1], tolower(data_vcf$ALT[which(data_vcf$POS == pos_name)]), tohama_full_genome[pos_numeric - 1])
        }
        if(tmp4%%3 == 0){
          tmp5 = paste0(tohama_full_genome[pos_numeric + 2], tohama_full_genome[pos_numeric + 1], tohama_full_genome[pos_numeric])
          tmp6 = paste0(tohama_full_genome[pos_numeric + 2], tohama_full_genome[pos_numeric + 1], tolower(data_vcf$ALT[which(data_vcf$POS == pos_name)]))
        }
        tmp5 = paste0(seqinr::comp(seqinr::s2c(tmp5)), collapse = '')
        tmp6 = paste0(seqinr::comp(seqinr::s2c(tmp6)), collapse = '')
      }
      if(data_genes$strand[a] == '-'){
        tmp4 = pos_numeric - data_genes$end[a] + 1
        if(tmp4%%3 == 1){
          tmp5 = paste0(tohama_full_genome[pos_numeric], tohama_full_genome[pos_numeric + 1], tohama_full_genome[pos_numeric + 2])
          tmp6 = paste0(tolower(data_vcf$ALT[which(data_vcf$POS == pos_name)]), tohama_full_genome[pos_numeric + 1], tohama_full_genome[pos_numeric + 2])
        }
        if(tmp4%%3 == 2){
          tmp5 = paste0(tohama_full_genome[pos_numeric - 1], tohama_full_genome[pos_numeric], tohama_full_genome[pos_numeric + 1])
          tmp6 = paste0(tohama_full_genome[pos_numeric - 1], tolower(data_vcf$ALT[which(data_vcf$POS == pos_name)]), tohama_full_genome[pos_numeric + 1])
        }
        if(tmp4%%3 == 0){
          tmp5 = paste0(tohama_full_genome[pos_numeric - 2], tohama_full_genome[pos_numeric - 1], tohama_full_genome[pos_numeric])
          tmp6 = paste0(tohama_full_genome[pos_numeric - 2], tohama_full_genome[pos_numeric - 1], tolower(data_vcf$ALT[which(data_vcf$POS == pos_name)]))
        }
      }
      possible_snps_genes[j] = tmp3
      if(is.na(tmp3)){
        possible_snps_genes[j] = tmp
      }
      possible_snps_pos_in_gene[j] = tmp4
      possible_snps_codon_ref[j] = tmp5
      possible_snps_codon_snp[j] = tmp6
      possible_snps_product[j] = data_genes$product[a]

      possible_snps_gene_ID_old[j] = tmp2
      idx = which(data_function_genes$Locus_id == tmp2)
      if(length(idx) >0){
        possible_snps_gene_function[j] = data_function_genes$Functional.category[idx]
        possible_snps_gene_prod_localisation[j] = data_function_genes$Subcellular.localization[idx]
      }
    }
    if(length(a)==0){
      if(is.na(pos_numeric) == F){
        s_min = which.min(abs(data_genes$end - pos_numeric)) ## Closest starting gene
        e_min = which.min(abs(data_genes$start - pos_numeric)) ## Closest end gene

        possible_snps_close_downstream_gene[j] = data_genes$gene[s_min]
        if(is.na(data_genes$gene[s_min])){
          possible_snps_close_downstream_gene[j] = data_genes$old_locus_tag[s_min]
        }
        possible_snps_close_downstream_gene_distance[j] = abs(data_genes$end - pos_numeric)[s_min]
        possible_snps_close_downstream_gene_strand[j] = data_genes$strand[s_min]

        possible_snps_close_upstream_gene[j] = data_genes$gene[e_min]
        if(is.na(data_genes$gene[e_min])){
          possible_snps_close_upstream_gene[j] = data_genes$old_locus_tag[e_min]
        }
        possible_snps_close_upstream_gene_distance[j] = abs(data_genes$start - pos_numeric)[e_min]
        possible_snps_close_upstream_gene_strand[j] = data_genes$strand[e_min]
      }

      tmp = NA
      tmp2 = NA
      tmp3 = NA
      tmp4 = NA

      if(possible_snps[j] == 'ptxP' | possible_snps[j] == 'fim3' | possible_snps[j] == 'prn'){
        possible_snps_genes[j] = possible_snps[j]
      }
    }
  }

  ## Save the reconstruction
  names_possibles_snps = data.frame('SNP' = possible_snps,
                                    'gene_name' = possible_snps_genes,
                                    'gene_id' = possible_snps_gene_ID_old,
                                    'gene_function' = possible_snps_gene_function,
                                    'gene_prod_localisation' = possible_snps_gene_prod_localisation,
                                    'product' = possible_snps_product,
                                    'pos_in_gene' = possible_snps_pos_in_gene,
                                    'codon_ref' = possible_snps_codon_ref,
                                    'codon_snp' = possible_snps_codon_snp)
  names_possibles_snps$AA_ref = translation_table[,2][match(names_possibles_snps$codon_ref, tolower(translation_table[,1]))]
  names_possibles_snps$AA_snp = translation_table[,2][match(names_possibles_snps$codon_snp, tolower(translation_table[,1]))]
  names_possibles_snps$name_snp = paste0(names_possibles_snps$AA_ref, floor((names_possibles_snps$pos_in_gene-1)/3)+1, names_possibles_snps$AA_snp)
  names_possibles_snps$syn = unlist(lapply(1:nrow(names_possibles_snps), function(x) names_possibles_snps$AA_ref[x] == names_possibles_snps$AA_snp[x] ))
  names_possibles_snps$close_upstream_gene = possible_snps_close_upstream_gene
  names_possibles_snps$close_upstream_gene_strand = possible_snps_close_upstream_gene_strand
  names_possibles_snps$close_upstream_gene_distance = possible_snps_close_upstream_gene_distance
  names_possibles_snps$close_downstream_gene = possible_snps_close_downstream_gene
  names_possibles_snps$close_downstream_gene_strand = possible_snps_close_downstream_gene_strand
  names_possibles_snps$close_downstream_gene_distance = possible_snps_close_downstream_gene_distance

  return(names_possibles_snps)
}

reconstruction_France$names_possibles_snps = annotate_possible_snps_tohama(reconstruction_France$possible_snps, data_vcf = data_vcf_France)

## Name promoter snps
promoter_snps = which(reconstruction_France$names_possibles_snps$close_downstream_gene_distance <= 100)
reconstruction_France$names_possibles_snps$syn[promoter_snps] = 'promoter'
reconstruction_France$names_possibles_snps$gene_name[promoter_snps] =  paste0('P', reconstruction_France$names_possibles_snps$close_downstream_gene[promoter_snps])
for(i in 1:length(promoter_snps)){
  idx = which(reconstruction_France$names_possibles_snps$gene_name == reconstruction_France$names_possibles_snps$close_downstream_gene[promoter_snps[i]])
  if(length(idx) > 0){
    reconstruction_France$names_possibles_snps$gene_function[promoter_snps[i]] = reconstruction_France$names_possibles_snps$gene_function[idx[1]]
  }
}

promoter_snps = which(reconstruction_France$names_possibles_snps$close_upstream_gene_distance <= 100)
reconstruction_France$names_possibles_snps$syn[promoter_snps] = 'promoter'
reconstruction_France$names_possibles_snps$gene_name[promoter_snps] = paste0('P', reconstruction_France$names_possibles_snps$close_upstream_gene[promoter_snps])
for(i in 1:length(promoter_snps)){
  idx = which(reconstruction_France$names_possibles_snps$gene_name == reconstruction_France$names_possibles_snps$close_upstream_gene[promoter_snps[i]])
  if(length(idx) > 0){
    reconstruction_France$names_possibles_snps$gene_function[promoter_snps[i]] = reconstruction_France$names_possibles_snps$gene_function[idx[1]]
  }
}


## Save output
# saveRDS(reconstruction_France, file = 'France_reconstruction_annotated_11032024.rds')
########################################################################################################################################












########################################################################################################################################
## PART 2: find linaege defining snps
#######################################################################################################################################
## Read output with all SNPs >0.1% and < 99.9% - this has been computed with the Part 1, and min_prop = 0.01 and max_prop = 0.99
reconstruction_France = readRDS(file = '3_BordetellaPertussis/2_analysis_index/3_snp_association/France_reconstruction_annotated_001_099.rds')

#####################################################################################
## Associate SNPs to groups
#####################################################################################
association_scores_per_group = function(dataset_with_nodes, dataset_with_inferred_reconstruction, tree, 
                                        possible_snps, names_possibles_snps,
                                        upstream_window, downstream_window){
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
      if(!is.na(dataset_with_inferred_reconstruction[tree$Nnode+2, which(colnames(dataset_with_inferred_reconstruction) == possible_snps[i])])){
        branches_snp = branches
        branches_snp[,1] = dataset_with_inferred_reconstruction[branches[,1],
                                                                which(colnames(dataset_with_inferred_reconstruction) == possible_snps[i])]
        branches_snp[,2] = dataset_with_inferred_reconstruction[branches[,2],
                                                                which(colnames(dataset_with_inferred_reconstruction) == possible_snps[i])]
        
        # ancestral_state = branches_snp[which.min(branches[,1]),1]
        
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
        
        ancestral_state = round(ancestral_state)
        branches_snp = round(branches_snp)
        branches_snp_bin = (branches_snp!=ancestral_state)
        
        tmp = which(branches_group[,1] == 1 & branches_group[,2] == 1)
        Px = as.numeric(branches_group[tmp,])
        Sx = as.numeric(branches_snp_bin[tmp,])
        score3 <- (sum((1 - Px)*(1 - Sx), na.rm=TRUE) + sum(Px*Sx, na.rm=TRUE)) / length(Px)
        
        scores[[j]][i] = score3
        
        t = table(branches_snp)
        t = t[which(names(t) != ancestral_state)]
        t = t[which.max(t)]
      }else{
        scores[[j]][i] = NA
      }
    }
    scores[[j]] = scores[[j]]#- median(scores[[j]])
    names(scores[[j]]) = possible_snps
  }
  return(scores)
}
## Time window
upstream_window = 10 # years (going to the root)
downstream_window = 100 # years (considering all the group)

## France
## Association scores
scores_France = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_France,
                                             dataset_with_inferred_reconstruction = reconstruction_France$dataset_with_inferred_resonstruction, 
                                              tree = tree_France, 
                                              possible_snps = reconstruction_France$possible_snps, 
                                              names_possibles_snps = reconstruction_France$names_possibles_snps$SNP,
                                              upstream_window, downstream_window)
## Find significant snps
edge_lineage_tree = split_France$lineage_tree$edge
edge_lineage_tree[,1] = split_France$tip_and_nodes_groups[match(edge_lineage_tree[,1],names(split_France$tip_and_nodes_groups))]
edge_lineage_tree[,2] = split_France$tip_and_nodes_groups[match(edge_lineage_tree[,2],names(split_France$tip_and_nodes_groups))]
edge_lineage_tree_snps = as.list(edge_lineage_tree)
scores_France_sig = scores_France
scores_France_sig_names = scores_France
sig_threshold = 0.8
for(i in 1:length(scores_France)){
  tmp = as.numeric(edge_lineage_tree[which(edge_lineage_tree[,2] == i),])
  if(tmp[1] == length(scores_France)+1) {
    ans = NULL
  }else{
    ans = names(which(scores_France[[tmp[1]]] > sig_threshold))
  }
  des = names(which(scores_France[[tmp[2]]] > sig_threshold))
  
  res = des[which(is.na(match(des, ans)))]
  scores_France_sig[[i]] = res
  
  scores_France_sig_names[[i]] = reconstruction_France$names_possibles_snps[match(res, reconstruction_France$names_possibles_snps$SNP),]
}
#####################################################################################

#####################################################################################
## First appearance groups
#####################################################################################
first_appearence_group_France = rep(NA, length(unique(dataset_with_nodes_France$groups)))
for(i in 1:length(first_appearence_group_France)){
  tmp = min(dataset_with_nodes_France$time[which(dataset_with_nodes_France$groups == i)])
  first_appearence_group_France[i] = tmp
}
#####################################################################################

#####################################################################################
## Remove lineage-defining SNPs before 1950
#####################################################################################
for(i in 1:length(scores_France_sig)){
  if(first_appearence_group_France[i] < 1960){
    scores_France_sig[[i]] = which(is.na(match(NULL, 0)))
    scores_France_sig_names[[i]] = reconstruction_France$names_possibles_snps[which(is.na(match(NULL, 0))),]
  }
}
#####################################################################################

#####################################################################################
## List significant snps
#####################################################################################
for(i in 1:length(scores_France_sig_names)){
  scores_France_sig_names[[i]]$group = rep(i, nrow(scores_France_sig_names[[i]]))
}
scores_France_sig_names_df = do.call("rbind",scores_France_sig_names)
#####################################################################################

#####################################################################################
## Create meaningful lineage trees (Figure 4c)
#####################################################################################
create_timed_lineage_tree = function(condensed_tree, node_group, dataset, time_root = 1900, tip_date = NULL){
  if(length(which(is.na(node_group))) > 0){
    tip_to_remove = which(is.na(condensed_tree$tip.label))
    condensed_tree$tip.label = condensed_tree$tip.label[-tip_to_remove]
    condensed_tree$edge = condensed_tree$edge[-which(condensed_tree$edge[,2] == tip_to_remove),]
    condensed_tree$edge[which(condensed_tree$edge<= length(condensed_tree$tip.label), arr.ind = T)] = condensed_tree$edge[which(condensed_tree$edge<= length(condensed_tree$tip.label), arr.ind = T)]-1
    condensed_tree$edge[which(condensed_tree$edge> length(condensed_tree$tip.label), arr.ind = T)] = condensed_tree$edge[which(condensed_tree$edge> length(condensed_tree$tip.label), arr.ind = T)]-1
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
lineage_tree_France = create_timed_lineage_tree(split_France$lineage_tree, 
                                                split_France$tip_and_nodes_groups, 
                                                dataset_with_nodes_France)

groups = levels(factor(split_France$tip_and_nodes_groups))
sig_snp_nodes_France_lineages = rep(NA, length(lineage_tree_France$node_group))
for(i in 1:(length(groups)-1)){
  idx = which(as.character(lineage_tree_France$node_group) == groups[i])
  tmp = scores_France_sig_names[[i]]
  tmp = tmp[which(tmp$syn == F | tmp$syn == 'promoter'),]
  tmp$gene_name[which(is.na(tmp$gene_name))] = tmp$SNP[which(is.na(tmp$gene_name))]
  if(length(idx) == 1){
    sig_snp_nodes_France_lineages[idx[1]] = paste0(tmp$gene_name, collapse = '/')
  }else if(length(idx) == 2){
    sig_snp_nodes_France_lineages[idx[2]] = paste0(tmp$gene_name, collapse = '/')
    sig_snp_nodes_France_lineages[idx[1]] = ''
  }else if(length(idx) > 2){
    sig_snp_nodes_France_lineages[idx[2]] = paste0(tmp$gene_name, collapse = '/')
    sig_snp_nodes_France_lineages[idx[1]] = ''
    sig_snp_nodes_France_lineages[idx[3:length(idx)]] = ''
  }
}
sig_snp_nodes_France_lineages[which(is.na(sig_snp_nodes_France_lineages))] = ''

order_France = order(unique(as.numeric(lineage_tree_France$node_group)))
p_France = ggtree(lineage_tree_France$lineage_tree, aes(col = lineage_tree_France$node_group), layout="roundrect", size=0.5, mrsd="2020-01-01") +
  geom_nodepoint(size=1.5, alpha=1, aes(col = lineage_tree_France$node_group)) +
  geom_tippoint(size=1.5, alpha=1, aes(col = lineage_tree_France$node_group)) +
  geom_text(size=1.5, aes(label = sig_snp_nodes_France_lineages), hjust = 1, vjust = -2)+
  scale_color_manual(values=colors_groups_France, breaks = unique(as.numeric(lineage_tree_France$node_group))[order_France], na.value = 'white', name = 'Groups') +
  theme_tree2()
p_France   
#####################################################################################

#####################################################################################
## Plot association scores (Figure 4g)
#####################################################################################
plot_scores_snp_France = function(){
  par(oma = c(0,0,0,0), mar = c(2,2,0,0), mgp = c(0,0.1,0), family = 'Arial',
      cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)
  colors = rep('grey', length(reconstruction_France$names_possibles_snps$SNP))
  idx = match(reconstruction_France$names_possibles_snps$SNP, rownames(reconstruction_France$snp_prevalence))
  colors[which(reconstruction_France$snp_prevalence[idx,3] > 0.001 & reconstruction_France$snp_prevalence[idx,3] < 1-0.001)] = 'black'
  plot(NULL, 
       xlim = c(1, max(as.numeric(lapply(reconstruction_France$names_possibles_snps$SNP, function(x)str_split(x, '_')[[1]][1])), na.rm = T)),
       ylim = c(0,1), bty = 'n', yaxt = 'n')
  axis(2, las = 2)
  tmp = apply(do.call(rbind, scores_France), MAR = 2, max)
  idx = which(reconstruction_France$names_possibles_snps$syn == F | 
                reconstruction_France$names_possibles_snps$syn == 'promoter')
  pos = as.numeric(lapply(reconstruction_France$names_possibles_snps$SNP[idx], function(x)str_split(x, '_')[[1]][1]))
  tmp = tmp[idx]
  points(pos, tmp, pch = 16, cex = 0.5, col = colors)
  abline(h = 0.8, lty =2, col = 'firebrick')
}
ggdraw(plot_scores_snp_France)
########################################################################################################################################

#####################################################################################
## Group by function (Figure 4k) 
#####################################################################################
plot_hist_gene_functions_France = function(){
  par(oma = c(0,0,0,0), mar = c(2,2,0,0), mgp = c(0,0.1,0), family = 'Arial',
      cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)
  
  t_function_all_snps_France = table(reconstruction_France$names_possibles_snps$gene_function[which(reconstruction_France$names_possibles_snps$syn == F |
                                                                                                      reconstruction_France$names_possibles_snps$syn == 'promoter')])
  t_function_all_snps_France = t_function_all_snps_France[which(names(t_function_all_snps_France) != 'pseudogenes' & names(t_function_all_snps_France) != 'pseudogenes (phase-variable)')]
  gene_cats = names(t_function_all_snps_France)
  
  scores_sig = do.call(rbind, scores_France_sig_names)
  scores_sig = scores_sig[which(scores_sig$syn == F | scores_sig$syn == 'promoter'),]

  t_snp_gene_function_France = table(factor(scores_sig$gene_function, levels = gene_cats))
  n = NULL
  for(i in 1:length(t_snp_gene_function_France)){
    n = c(n, t_function_all_snps_France[which(names(t_function_all_snps_France) == names(t_snp_gene_function_France[i]))])
  }

  t_snp_gene_function_France = table(factor(scores_sig$gene_function, levels = gene_cats))
  
  for(i in 1:length(t_snp_gene_function_France)){
    t_snp_gene_function_France[i] = t_snp_gene_function_France[i]/t_function_all_snps_France[which(names(t_function_all_snps_France) == names(t_snp_gene_function_France[i]))]
  }
  
  t_snp_gene_function_France = t_snp_gene_function_France[order(t_snp_gene_function_France, decreasing = T)]
  t_snp_gene_function_France = t_snp_gene_function_France[which(t_snp_gene_function_France > 0)]
  
  plot(x = 1:(length(t_snp_gene_function_France)+1),
       y = c(t_snp_gene_function_France,0), 
       ylim = c(0,max(t_snp_gene_function_France)), col = c('firebrick', rep('black', length(t_snp_gene_function_France)-1), 'grey'),
       type = 'h', lwd = 5, bty = 'n', xaxt = 'n', yaxt = 'n',
       ylab = 'Proportion', xlab = '')
  axis(2, las = 2)
  axis(1, at = seq(1, length(t_snp_gene_function_France)+1, 1), labels = rep('', length(t_snp_gene_function_France)+1))
  text(x = 1:(length(t_snp_gene_function_France)+1),
       y = par("usr")[3] - 0.001,
       labels = c(names(t_snp_gene_function_France), 'all others'),
       xpd = NA,
       adj = 1,
       cex = 0.3,
       ## Rotate the labels by 35 degrees.
       srt = 35)
}
ggdraw(plot_hist_gene_functions_France)
#####################################################################################

#####################################################################################
## Save multipanel
#####################################################################################
panel_France = plot_grid(p_France, ggdraw(plot_scores_snp_France), ggdraw(plot_hist_gene_functions_France),
                  rel_widths = c(1, 1.5, 1.25), labels = '', ncol = 3)
ggsave(filename = '3_BordetellaPertussis/2_analysis_index/3_snp_association/Panel_association_Pertussis_France.pdf', 
       plot = panel_France, device = 'pdf', scale = 1,
       width = 22, height = 5, units = 'cm')
########################################################################################################################################


