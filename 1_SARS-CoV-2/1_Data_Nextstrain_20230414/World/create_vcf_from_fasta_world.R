## Create VCF like file for GWAS
library(seqinr)
setwd('/Users/noemielefrancq/Documents/THD/THD_SARS-CoV-2/1_Data_Nextstrain_20230414/World/')

## Input fasta data
fasta_data = read.fasta('nextstrain_ncov_gisaid_global_all-time_timetree_dedup_wuhan.fasta', forceDNAtolower = T, whole.header = T)
data_seq = getSequence(fasta_data)
data_name_seq = names(fasta_data)

## Transform into matrix of characters
matrix_data = matrix(0, length(data_name_seq), length(data_seq[[1]]))
for (i in 1:length(data_seq)){
  matrix_data[i,] =  data_seq[[i]]
}
rownames(matrix_data) = data_name_seq
matrix_data_nuc = matrix_data
remove(fasta_data)

## Split this matrix into ORFs, and translate to AA
data_orfs = read.csv('../../Reference/Position_ORFs_nextstrain.csv')
list_matrices_ORFs = list_matrices_ORFs_codon = list()
for(i in 1:nrow(data_orfs)){
  # list_matrices_ORFs[[i]] = t(apply(matrix_data[,data_orfs$start[i]:data_orfs$end[i]], MARGIN = 1, function(x)translate(x)))
  
  ## AA sequences
  list_matrices_ORFs[[i]] = t(apply(matrix_data[,data_orfs$start[i]:data_orfs$end[i]], MARGIN = 1, function(x)translate(x)))
  ## Matrix of codons
  list_matrices_ORFs_codon[[i]] = matrix(NA, nrow = nrow(matrix_data),
                                         ncol = ncol(matrix_data[,data_orfs$start[i]:data_orfs$end[i]])/3)
  for(j in 1:(ncol(list_matrices_ORFs_codon[[i]]))){
    list_matrices_ORFs_codon[[i]][,j] = apply(matrix_data[,data_orfs$start[i]:data_orfs$end[i]][,(1:3)+(j-1)*(3)], 
                                              MARGIN = 1, function(x)paste0(x, collapse = ''))
  }
  ## Replace gaps and IUPAC code by N
  list_matrices_ORFs_codon[[i]] = t(apply(list_matrices_ORFs_codon[[i]], 
                                          MARGIN = 1, 
                                          FUN = function(x)str_replace_all(string = x, pattern = 'r|y|s|w|k|m|b|d|h|v|n|-', replacement = 'n')))
  ## Translate again to check
  for(j in 1:ncol(list_matrices_ORFs_codon[[i]])){
    for(k in 1:length(list_matrices_ORFs_codon[[i]][j,])) {
      list_matrices_ORFs_codon[[i]][j,k] = translate(str_split(list_matrices_ORFs_codon[[i]][j,k], pattern = '')[[1]])
    }
  }
}

## Function to transform matrix into binary matrix
bin <- function(x){
  x=as.numeric(factor(as.character(x),levels=c("A","C","D","E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", 
                                               "V", "W", "X", "Y", "*")))
  return(x)
}
bin_nuc <- function(x){
  x=as.numeric(factor(as.character(x),levels=c("a","g","c","t", "n", "-", "b", "d", "h", "v", "r", "y", "s", "w", "k", "m")))
  return(x)
}

###########################################################3
## For each ORF, find SNPs and wirte them in a vcf
for(orf in 1:length(list_matrices_ORFs)){
  matrix_data = list_matrices_ORFs[[orf]]
  
  ## X: unsure site, replace by ref
  for(i in 1:ncol(matrix_data)){
    matrix_data[which(matrix_data[,i] == "X"),i] = matrix_data[1,i]
  }
  
  matrix_data_bin = apply(matrix_data, MARGIN = 2, FUN=bin)
  
  ## Find SNPs
  matrix_data_bin_snp = t(apply(matrix_data_bin, MARGIN = 1, FUN=function(x)x-matrix_data_bin[1,]))
  matrix_data_bin_snp = apply(matrix_data_bin_snp, MARGIN = 2, FUN=function(x)as.numeric(factor(as.numeric(factor(x, levels = c(x[1], -30:(x[1]-1), (x[1]+1):(30))))))-1)
  matrix_data_bin_snp[1:10, 1:10]
  
  ## Subset matrix, where there are snps
  overall_snps = colSums(matrix_data_bin_snp)
  a = which(overall_snps > 0)
  pos_snp = (1:ncol(matrix_data_bin_snp))[a]
  matrix_data_bin_snp_only = matrix_data_bin_snp[,a]
  matrix_data_bin_snp_only[1:10, 1:10]
  dim(matrix_data_bin_snp_only)
  
  ## Find positions where there are 1+ SNPs
  matrix_data_bin_snp_only_singles = NULL
  pos_snp_singles = NULL
  pos_snp_singles_names = NULL
  for(i in 1:ncol(matrix_data_bin_snp_only)){
    if(i %% 100 == 0){
      print(paste0(i, ' / ', ncol(matrix_data_bin_snp_only)))
    }
    
    t = table(matrix_data_bin_snp_only[,i])
    if(length(t) <= 2){
      if(!is.null(matrix_data_bin_snp_only_singles)) {
        matrix_data_bin_snp_only_singles = cbind(matrix_data_bin_snp_only_singles, matrix_data_bin_snp_only[,i])
      } 
      if(is.null(matrix_data_bin_snp_only_singles)) {
        matrix_data_bin_snp_only_singles = matrix_data_bin_snp_only[,i]
      }
      pos_snp_singles = c(pos_snp_singles, pos_snp[i])
      pos_snp_singles_names = c(pos_snp_singles_names, pos_snp[i])
    }
    if(length(t) > 2){
      # print(t)
      for(j in 2:length(t)){
        tmp = matrix_data_bin_snp_only[,i]
        tmp[which(tmp != (j-1))] = 0
        tmp[which(tmp > 0)] = 1
        matrix_data_bin_snp_only_singles = cbind(matrix_data_bin_snp_only_singles, tmp)
      }
      pos_snp_singles_names = c(pos_snp_singles_names, paste0(rep(pos_snp[i], length(t)-1), '_', LETTERS[1:length(t)-1]))
      pos_snp_singles = c(pos_snp_singles, rep(pos_snp[i], length(t)-1))
    }
  }
  colnames(matrix_data_bin_snp_only_singles) = pos_snp_singles
  
  ## Write genotype vcf
  data_vcf = data.frame('#CHROM' = rep(paste0(data_orfs$gene[orf]), dim(matrix_data_bin_snp_only_singles)[2]),
                        'POS' = pos_snp_singles_names,
                        'ID' = rep('.', dim(matrix_data_bin_snp_only_singles)[2]),
                        'REF' = rep('NA', dim(matrix_data_bin_snp_only_singles)[2]),
                        'ALT' = rep('NA', dim(matrix_data_bin_snp_only_singles)[2]),
                        'QUAL'  = rep('5000', dim(matrix_data_bin_snp_only_singles)[2]),
                        'FILTER'  = rep('PASS', dim(matrix_data_bin_snp_only_singles)[2]),
                        'INFO'  = rep('.', dim(matrix_data_bin_snp_only_singles)[2]),
                        'FORMAT' = rep('GT', dim(matrix_data_bin_snp_only_singles)[2]))
  for(i in 1:dim(matrix_data_bin_snp_only_singles)[1]){
    data_vcf = cbind(data_vcf, c(matrix_data_bin_snp_only_singles[i,]))
  }
  for(i in 1:dim(matrix_data_bin_snp_only_singles)[2]){
    data_vcf$REF[i] = toupper(matrix_data[1,pos_snp_singles[i]])
    t = table(matrix_data[which(matrix_data_bin_snp_only_singles[,i] != 0),pos_snp_singles[i]])
    snp = toupper(names(t))
    if(length(t)>1){
      print(snp)
      print(i)
    }
    data_vcf$ALT[i] = snp
  }
  colnames(data_vcf) = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', data_name_seq)
  
  ## Write output
  write.table(data_vcf, file = paste0('nextstrain_ncov_gisaid_global_all-time_timetree_', data_orfs$gene[orf], '.vcf'), quote = F, na = '.', append = F, sep = '\t', row.names = F, 
              col.names = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', data_name_seq))
}
###########################################################

###########################################################
## Full genome nucleotide changes
matrix_data = matrix_data_nuc

## X: unsure site, replace by ref
for(i in 1:ncol(matrix_data)){
  matrix_data[which(matrix_data[,i] == "-" | matrix_data[,i] == "n"),i] = matrix_data[1,i]
}

matrix_data_bin = apply(matrix_data, MARGIN = 2, FUN=bin_nuc)

## Find SNPs
matrix_data_bin_snp = t(apply(matrix_data_bin, MARGIN = 1, FUN=function(x)x-matrix_data_bin[1,]))
matrix_data_bin_snp = apply(matrix_data_bin_snp, MARGIN = 2, FUN=function(x)as.numeric(factor(as.numeric(factor(x, levels = c(x[1], -30:(x[1]-1), (x[1]+1):(30))))))-1)
matrix_data_bin_snp[1:10, 1:10]

## Subset matrix, where there are snps
overall_snps = colSums(matrix_data_bin_snp)
a = which(overall_snps > 0)
pos_snp = (1:ncol(matrix_data_bin_snp))[a]
matrix_data_bin_snp_only = matrix_data_bin_snp[,a]
matrix_data_bin_snp_only[1:10, 1:10]
dim(matrix_data_bin_snp_only)

## Find positions where there are 1+ SNPs
matrix_data_bin_snp_only_singles = NULL
pos_snp_singles = NULL
pos_snp_singles_names = NULL
for(i in 1:ncol(matrix_data_bin_snp_only)){
  if(i %% 100 == 0){
    print(paste0(i, ' / ', ncol(matrix_data_bin_snp_only)))
  }
  
  t = table(matrix_data_bin_snp_only[,i])
  if(length(t) <= 2){
    if(!is.null(matrix_data_bin_snp_only_singles)) {
      matrix_data_bin_snp_only_singles = cbind(matrix_data_bin_snp_only_singles, matrix_data_bin_snp_only[,i])
    } 
    if(is.null(matrix_data_bin_snp_only_singles)) {
      matrix_data_bin_snp_only_singles = matrix_data_bin_snp_only[,i]
    }
    pos_snp_singles = c(pos_snp_singles, pos_snp[i])
    pos_snp_singles_names = c(pos_snp_singles_names, pos_snp[i])
  }
  if(length(t) > 2){
    # print(t)
    for(j in 2:length(t)){
      tmp = matrix_data_bin_snp_only[,i]
      tmp[which(tmp != (j-1))] = 0
      tmp[which(tmp > 0)] = 1
      matrix_data_bin_snp_only_singles = cbind(matrix_data_bin_snp_only_singles, tmp)
    }
    pos_snp_singles_names = c(pos_snp_singles_names, paste0(rep(pos_snp[i], length(t)-1), '_', LETTERS[1:length(t)-1]))
    pos_snp_singles = c(pos_snp_singles, rep(pos_snp[i], length(t)-1))
  }
}
colnames(matrix_data_bin_snp_only_singles) = pos_snp_singles

## Write genotype vcf
data_vcf = data.frame('#CHROM' = rep(paste0('nuc'), dim(matrix_data_bin_snp_only_singles)[2]),
                      'POS' = pos_snp_singles_names,
                      'ID' = rep('.', dim(matrix_data_bin_snp_only_singles)[2]),
                      'REF' = rep('NA', dim(matrix_data_bin_snp_only_singles)[2]),
                      'ALT' = rep('NA', dim(matrix_data_bin_snp_only_singles)[2]),
                      'QUAL'  = rep('5000', dim(matrix_data_bin_snp_only_singles)[2]),
                      'FILTER'  = rep('PASS', dim(matrix_data_bin_snp_only_singles)[2]),
                      'INFO'  = rep('.', dim(matrix_data_bin_snp_only_singles)[2]),
                      'FORMAT' = rep('GT', dim(matrix_data_bin_snp_only_singles)[2]))
for(i in 1:dim(matrix_data_bin_snp_only_singles)[1]){
  data_vcf = cbind(data_vcf, c(matrix_data_bin_snp_only_singles[i,]))
}
for(i in 1:dim(matrix_data_bin_snp_only_singles)[2]){
  data_vcf$REF[i] = toupper(matrix_data[1,pos_snp_singles[i]])
  t = table(matrix_data[which(matrix_data_bin_snp_only_singles[,i] != 0),pos_snp_singles[i]])
  snp = toupper(names(t))
  if(length(t)>1){
    print(snp)
    print(i)
  }
  data_vcf$ALT[i] = snp
}
colnames(data_vcf) = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', data_name_seq)

## Write output
write.table(data_vcf, file = paste0('nextstrain_ncov_gisaid_global_all-time_timetree_', 'full', '.vcf'), quote = F, na = '.', append = F, sep = '\t', row.names = F, 
            col.names = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', data_name_seq))
###########################################################











