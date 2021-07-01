#!/bin/bash
#SBATCH --nodes=1
#SBATCH --array=0-58:1
#SBATCH --time=120:00:00
#SBATCH --job-name=ANGSD_SNP_Calling
#SBATCH --output=ANGSD_SNP_Calling_%A_%a.log
#SBATCH --error=ANGSD_SNP_Calling_%A_%a.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@zoo.ox.ac.uk

#########################################################################################################
# CONDUCT GENOTYPE LIKELIHOOD ESTIMATION WITH ANGSD v.0.921 AND IMPUTE AND PHASE GENOTYPES USING BEAGLE4
# Will output the following files per scaffold:
# 1. VCF with genotype probability (GP) and genotype likelihood (GL) fields
# 2. MAFs
# 3. Genotype likelihoods (GLs) in beagle format
# 4. VCF with imputed and phased genotypes

# A. Sendell-Price, June 2021
#########################################################################################################

#Load required modules
module load angsd/0.921
module load java/1.8.0
module load vcftools

# Set path to reference assembly and list of bam files (bam.list)
# Note: bam files need to be indexed (using samtools index) 
REF=/data/zool-zost/Ref_Genome/GCA_001281735.1_ASM128173v1_genomic.fna
BAMs=bam.list

# Create direcories for output (ignored (with warning) if already exists)
mkdir MAFs VCFs BEAGLEs BEAGLE_LOGs
mkdir VCFs/RAW
mkdir VCFs/IMPUTED
mkdir VCFs/PHASED

# Use slurm array task ID to get list of scaffold names
SCAFFOLDS=scaffold_lists/scaffold.list.$(printf '%02d' $SLURM_ARRAY_TASK_ID)
END=$(cat $SCAFFOLDS | wc -l)

# For each line within scaffold list do the following
for i in $(eval echo "{1..$END}")
do

# Extact scaffold ID
REGION_2_LOOKUP=$(cat $SCAFFOLDS | head -n $i | tail -n 1)

# Estimate GLs and GPs using ANGSD
# Note: No attempt has been made to call an individuals genotype, instead these will be imputed later
angsd -b $BAMs -ref $REF \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -baq 1 -minMapQ 20 -minQ 20 \
-GL 1 -doMajorMinor 1 -doMaf 1 -doPost 2 -doVcf 1 -doGlf 2 -minMaf 0.05 -SNP_pval 1e-6 -minInd 194 \
-r $REGION_2_LOOKUP -out ZLat.v1_${REGION_2_LOOKUP}

# Explanation of above settings:
# ==============================
# -uniqueOnly = only use uniquely mapped reads (ingnore reads with multiple hits)
# -remove_bads = same as the samtools flags -x which removes read with a flag above 255 (not primary, failure and duplicate reads)
# -only_proper_pairs = include only pairs of reads with both mates (forward and reverse) mapped correctly
# -baq = perform BAQ computation (1: same as in SAMtools)
# -minMapQ = minimum mapQ quality
# -minQ = minimum base quality score
# -GL = calculate genotype likelihoods (1: using SAMtools model )
# -doMajorMinor = infer major and minor alleles (1: from GLs)
# -doPost = calculate posterior prob (1: Using frequency as prior)
# -doVcf = output a VCF file (1: yes)
# -doGlf = ouput genotype likelihoods file (4: beagle likelihood format)
# -minMaf = minumum minor allele frequency tolerated
# -SNP_pval = significance threshold for determining true polymorphism
# -minInd = only output sites with information for atleast [int] samples (half is sensible)
#           Although we will impute missing data, this means we will only impute sites where we
#			have a good idea what is going on.
# ==============================

# Use beagle4.1 to impute and phase genotypes
java -Xmx15000m -jar /data/zool-zost/cont7348/BIN/beagle.27Jan18.7e1.jar \
gl=ZLat.v1_${REGION_2_LOOKUP}.vcf.gz out=ZLat.v1_${REGION_2_LOOKUP}.imputed

java -Xmx15000m -jar /data/zool-zost/cont7348/BIN/beagle.27Jan18.7e1.jar \
gt=ZLat.v1_${REGION_2_LOOKUP}.imputed.vcf.gz out=ZLat.v1_${REGION_2_LOOKUP}.imputed.phased

# Perform some tidying up
mv ZLat.v1_${REGION_2_LOOKUP}.mafs.gz MAFs/
mv ZLat.v1_${REGION_2_LOOKUP}.vcf.gz VCFs/RAW/
mv ZLat.v1_${REGION_2_LOOKUP}.beagle.gz BEAGLEs/
mv ZLat.v1_${REGION_2_LOOKUP}.imputed.vcf.gz VCFs/IMPUTED/
mv ZLat.v1_${REGION_2_LOOKUP}.imputed.phased.vcf.gz VCFs/IMPUTED/
mv ZLat.v1_${REGION_2_LOOKUP}.imputed*.log BEAGLE_LOGs/
rm ZLat.v1_${REGION_2_LOOKUP}.arg

done
