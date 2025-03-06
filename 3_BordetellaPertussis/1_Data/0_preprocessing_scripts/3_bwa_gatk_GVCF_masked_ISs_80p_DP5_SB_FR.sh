#!/bin/sh

#########################################################
## BWA mapping + Variant calling with GATK
## from Nadia's script
## Adapted for Cambridge cluster
### LATEST ### 06-07-2021

## Modules require, but loaded before
#module load bwa-0.7.17-gcc-5.4.0-42mry2g
#module load samtools/1.9

#################################################### Mapping with BWA
#Create pe.sam files out of trimmed fastq
## Decompress trimmed fastqs
pigz -p 4 -d --keep fastq_trimmed_dedup/$( basename ${file} _1_trim_dedup.fastq.gz)_1_trim_dedup.fastq.gz ;
pigz -p 4 -d --keep fastq_trimmed_dedup/$( basename ${file} _1_trim_dedup.fastq.gz)_2_trim_dedup.fastq.gz ;

## Move files to current directory
mv fastq_trimmed_dedup/$( basename ${file} _1_trim_dedup.fastq.gz)_1_trim_dedup.fastq fastq_trimmed_dedup/$( basename ${file} _1_trim_dedup.fastq.gz)_2_trim_dedup.fastq ./

## Run bwa
R1=$( basename ${file} _1_trim_dedup.fastq.gz)""_1_trim_dedup.fastq
R2=$( basename ${file} _1_trim_dedup.fastq.gz)""_2_trim_dedup.fastq
/home/ncmjl2/softwares/bwa/bwa mem -t 4 /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta $R1 $R2 > $( basename ${file} _1_trim_dedup.fastq.gz)_trim.pe.sam

## Clean directory
rm $R1 ;
rm $R2 ;

##################################################### Variant calling
## Make tmp dir
#mkdir tmp_$( basename ${file} _1_trim_dedup.fastq.gz)
#export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp_$( basename ${file} _1_trim_dedup.fastq.gz)

## Sort sam file
/home/ncmjl2/softwares/gatk-4.2.0.0/gatk SortSam -I $( basename ${file} _1_trim_dedup.fastq.gz)_trim.pe.sam -O $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_sorted.bam -SO coordinate
rm $( basename ${file} _1_trim_dedup.fastq.gz)_trim.pe.sam

## Mark duplicates, add read group and create a bam index
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk MarkDuplicates -I $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_sorted.bam -O $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup.bam -M $( basename ${file} _1_trim_dedup.fastq.gz)""metrics.txt 
java -jar /home/ncmjl2/softwares/picard/build/libs/picard.jar AddOrReplaceReadGroups I=$( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup.bam O=$( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup2.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$( basename ${file} _1_trim_dedup.fastq.gz)
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk BuildBamIndex -I $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup2.bam

## Get first variant calling done
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk HaplotypeCaller -R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta -I $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup2.bam --base-quality-score-threshold 20 --minimum-mapping-quality 60 --mapping-quality-threshold-for-genotyping 60 -ploidy 1 -ERC GVCF --annotation StrandBiasBySample --annotation AlleleFraction -G AS_StandardAnnotation -O $( basename ${file} _1_trim_dedup.fastq.gz)_GATK.g.vcf --all-site-pls

## Compute coverage
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk CollectWgsMetrics --USE_FAST_ALGORITHM -CAP 500 --MINIMUM_BASE_QUALITY 20 --MINIMUM_MAPPING_QUALITY 60 --INCLUDE_BQ_HISTOGRAM -R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta -I $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup2.bam -O $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage.txt 

## Extract interesting line
# awk '$1 == 4086189' $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage.txt > $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage_line.txt

## Extract mean coverage, and sd
#cat $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage_line.txt | awk -F'\t' '{print $2}' > $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage_mean.txt
#cat $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage_line.txt | awk -F'\t' '{print $3}' > $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage_sd.txt
#mean=$( cat $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage_mean.txt)
#sd=$( cat $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage_sd.txt )
#threshold=$( echo "scale=1 ; $mean - 2.5*$sd" | bc )
#threshold=$( echo $threshold | awk '{print int($1+0.5)}' )
#f [ "$threshold" -lt "2" ]; then
#    threshold=2;
#i
#if [ "$threshold" -gt "10" ]; then
#    threshold=10;
#fi
threshold=5

## Clean
rm $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup.bam 
rm $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_sorted.bam
rm $( basename ${file} _1_trim_dedup.fastq.gz)""metrics.txt 
rm $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage_sd.txt
rm $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage_mean.txt
rm $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage_line.txt
rm -r tmp_$( basename ${file} _1_trim_dedup.fastq.gz)

## Move GVCF file
mv $( basename ${file} _1_trim_dedup.fastq.gz)""_GATK.g.vcf g_vcf_IS_masked/
mv $( basename ${file} _1_trim_dedup.fastq.gz)""_GATK.g.vcf.idx g_vcf_IS_masked/
mv $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage.txt g_vcf_IS_masked/


##################################################### Filter variants, to take out IS
## Genotype
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk --java-options "-Xmx10g -Xms10g" GenotypeGVCFs \
   -R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta \
   -V g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_GATK.g.vcf \
   --sample-ploidy 1 \
   -G AS_StandardAnnotation \
   -A AlleleFraction \
   -O g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_genotyped.vcf.gz
      
## Filter IS
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk --java-options "-Xmx10g -Xms10g" VariantFiltration \
   -V g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_genotyped.vcf.gz \
   -filter "QD < 2.0" --filter-name "QD2" \
   -filter "MQ < 20.0" --filter-name "MQ20" \
   -filter "DP < $threshold" --filter-name "DPthre" \
   -filter "SOR > 3.0"  --filter-name "SORthre" \
   --genotype-filter-expression "AF < 0.75" --genotype-filter-name "AFthre" \
   --exclude-intervals /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/IS_all_Tohama_pertussis_full_positions_transposase_inc.intervals \
   -O g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.vcf.gz
   
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk SelectVariants \
    -R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta \
    -V g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.vcf.gz \
    -O g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.vcf.gz \
    --exclude-filtered \
    --set-filtered-gt-to-nocall \
    --select-type-to-include SNP

##################################################### Compute mask, to put N is the reference
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk HaplotypeCaller -base-quality-score-threshold 20 --minimum-mapping-quality 60 --mapping-quality-threshold-for-genotyping 60 -R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta -I $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup2.bam -ploidy 1 -ERC BP_RESOLUTION  -G AS_StandardAnnotation  --annotation StrandBiasBySample -O g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_GATK_allsite.g.vcf

/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk GenotypeGVCFs \
    -V g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_GATK_allsite.g.vcf \
    -R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta \
    -G AS_StandardAnnotation \
    --all-sites \
    -O g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_allsite_raw.vcf
    
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk SelectVariants \
    -R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta \
    -V g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_allsite_raw.vcf \
    -O g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked1.vcf \
    -select "DP < $threshold || SOR > 3.0 || QD < 2.0 || MQ < 20.0" 
    
grep :0:0:0 g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_allsite_raw.vcf > g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked2.vcf
    
cat g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked1.vcf g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked2.vcf > g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked.vcf

rm g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked1.vcf 
rm g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked2.vcf

cat g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked.vcf |tail -n +36 |awk '{FS="\t";OFS="\t";print $1,$2-1,$2,$3, etc}' > g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked.bed

##################################################### Compute fasta file
#/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk FastaAlternateReferenceMaker  \
#--snp-mask g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked2.vcf.gz  \
#--snp-mask-priority TRUE \
#-R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta \
#-V g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.vcf.gz \
#-O g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_consensus.fasta

/home/ncmjl2/softwares/bcftools-1.12/bcftools consensus -f /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta --mask g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked.bed g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.vcf.gz > g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_consensus.fasta

/home/ncmjl2/softwares/bcftools-1.12/bcftools annotate -x INFO,^FORMAT/GT -O z g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.vcf.gz > g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.GT.vcf.gz
##for i in *_filtered.PASS.GT.vcf.gz; do tabix $i; done

## Clean
rm g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_genotyped.vcf.gz* 
rm g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.vcf.gz*
mv g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.vcf.gz* filtered_PASS_vcf/
mv g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked.* filtered_PASS_vcf/
rm g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_GATK_allsite*
rm g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_allsite*
mv g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_consensus* consensus_fasta/
mv g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.GT.vcf.gz filtered_PASS_no_GT_vcf/
mv $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup2.bam bam_dedup/
mv $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup2.bai bam_dedup/
