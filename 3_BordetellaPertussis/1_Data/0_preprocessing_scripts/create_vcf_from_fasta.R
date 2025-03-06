## Create VCF like file from fasta file
library(seqinr)
library(vcfR)
setwd('2_data_by_country/')

country = 'France'
name = 'FR_80p_dp5_SOR_QC_names_with_Tohama_recom_masked_snp'

## Input fasta data (alignment from snp-sites)
fasta_data = read.fasta(paste0(country, '/', name,'.fasta'))
vcf_data = read.vcfR(paste0(country, '/', name,'.vcf'))
vcf_data = vcf_data@fix

data_seq = getSequence(fasta_data)
data_name_seq = names(fasta_data)
position_snps = vcf_data[,2]
  
## Transform into matrix of characters
matrix_data = matrix(0, length(data_name_seq), length(data_seq[[1]]))
for (i in 1:length(data_seq)){
  matrix_data[i,] =  data_seq[[i]]
}
rownames(matrix_data) = data_name_seq
colnames(matrix_data) = position_snps
remove(fasta_data)
  
## Put Tohama in the first row: to compute SNP based on it
a = which(data_name_seq == "NC_002929.2")
matrix_data = rbind(matrix_data[a,],
                    matrix_data[-a,])
rownames(matrix_data) = c("NC_002929.2", rownames(matrix_data)[-1])

## Transform matrix into binary matrix
bin <- function(x){
  x=as.numeric(factor(as.character(x),levels=c("a","g","c","t", "n", "-", "b", "d", "h", "v", "r", "y", "s", "w", "k", "m")))
  return(x)
}
matrix_data_bin = apply(matrix_data, MARGIN = 2, FUN=bin)

## Find SNPs
matrix_data_bin_snp = t(apply(matrix_data_bin, MARGIN = 1, FUN=function(x)x-matrix_data_bin[1,]))
matrix_data_bin_snp = apply(matrix_data_bin_snp, MARGIN = 2, FUN=function(x)as.numeric(factor(as.numeric(factor(x, levels = c(x[1], -20:(x[1]-1), (x[1]+1):(20))))))-1)
matrix_data_bin_snp[1:10, 1:10]

## Subset matrix 
overall_snps = colSums(matrix_data_bin_snp)
a = which(overall_snps > 0)
pos_snp = (colnames(matrix_data_bin_snp))[a]
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
    # print(length(pos_snp_singles))
    # print(dim(matrix_data_bin_snp_only_singles))
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
data_vcf = data.frame('#CHROM' = rep('CHR1', dim(matrix_data_bin_snp_only_singles)[2]),
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
rownames(matrix_data) = str_replace_all(rownames(matrix_data), '\\|', '_')
colnames(data_vcf) = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', rownames(matrix_data))

## Write output
write.table(data_vcf, file = paste0(country, '/', name,'_withN.vcf'), quote = F, na = '.', append = F, sep = '\t', row.names = F, 
            col.names = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',  rownames(matrix_data)))
## REMEMBER TO ADD '##fileformat=VCFv4.2' HEADER!
# for i in seq_jc_uni1000_*; do cat line.txt $i > V_$i; done




