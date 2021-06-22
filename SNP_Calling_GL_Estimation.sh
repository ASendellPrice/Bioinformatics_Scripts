#!/bin/bash
#SBATCH --clusters=htc
#SBATCH --array=0-58:1
#SBATCH --partition=long
#SBATCH --time=120:00:00
#SBATCH --job-name=ANGSD_SNP_Calling
#SBATCH --output=ANGSD_SNP_Calling_%A_%a.log
#SBATCH --error=ANGSD_SNP_Calling_%A_%a.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@zoo.ox.ac.uk

######################################################################################
# CONDUCT GENOTYPE LIKELIHOOD ESTIMATION AND CALL SAMPLE GENOTYPES WITH ANGSD v.0.921
# Will output the following files per scaffold:
# 1. VCF
# 2. GENO
# 3. MAFs
# 4. Beagle (GLs)

# Requires angsd to be compiled within a conda environment (see notes at end)

# A. Sendell-Price, June 2021
######################################################################################

# Load Anaconda module
module load Anaconda2

# Load conda environment
source activate $DATA/angsd-env

# Set path to reference assembly and list of bam files (bam.list)
# Note: bam files need to be indexed (using samtools index) 
REF=/data/zool-zost/Ref_Genome/GCA_001281735.1_ASM128173v1_genomic.fna
BAMs=bam.list

# Create direcories for output (ignored (with warning) if already exists)
mkdir GENOs MAFs VCFs BEAGLEs

# Use slurm array task ID to get list of scaffold names
SCAFFOLDS=scaffold_lists/scaffold.list.$(printf '%02d' $SLURM_ARRAY_TASK_ID)
END=$(cat $SCAFFOLDS | wc -l)

# For each line within scaffold list do the following
for i in $(eval echo "{1..$END}")
do

# Extact scaffold ID
REGION_2_LOOKUP=$(cat $SCAFFOLDS | head -n $i | tail -n 1)

# Call SNPs (and estimate GLs) using ANGSD
angsd -b $BAMs -ref $REF \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -baq 1 -minMapQ 20 -minQ 20 -geno_minDepth 5 \
-doCounts 1 -GL 1 -doMajorMinor 1 -doMaf 1 -doPost 2 -doGeno 4 -doVcf 1 -doGlf 2 -minMaf 0.01 -SNP_pval 1e-6 \
-r $REGION_2_LOOKUP -out ZLat.v1_${REGION_2_LOOKUP}

# Explanation of above settings:
# ==============================
# -uniqueOnly = only use uniquely mapped reads (ingnore reads with multiple hits)
# -remove_bads = same as the samtools flags -x which removes read with a flag above 255 (not primary, failure and duplicate reads)
# -only_proper_pairs = include only pairs of reads with both mates (forward and reverse) mapped correctly
# -baq = perform BAQ computation (1: same as in SAMtools)
# -minMapQ = minimum mapQ quality
# -minQ = minimum base quality score
# -geno_minDepth = only call genotypes if read depth is at least [int] for that individual
# -doCounts = count occurances of As Ts Cs Gs
# -GL = calculate genotype likelihoods (1: using SAMtools model )
# -doMajorMinor = infer major and minor alleles (1: from GLs)
# -doPost = calculate posterior prob (1: Using frequency as prior)
# -doGeno = call genotypes (4: print the called genotype as AA, AC, AG)
# -doVcf = output a VCF file (1: yes)
# -doGlf = ouput genotype likelihoods (4: beagle likelihood format)
# -minMaf = minumum minor allele frequency tolerated
# -SNP_pval = significance threshold for determining true polymorphism

# Perform some tidying up
mv ZLat.v1_${REGION_2_LOOKUP}.geno.gz GENOs/
mv ZLat.v1_${REGION_2_LOOKUP}.mafs.gz MAFs/
mv ZLat.v1_${REGION_2_LOOKUP}.vcf.gz VCFs/
mv ZLat.v1_${REGION_2_LOOKUP}.beagle.gz BEAGLEs/
rm ZLat.v1_${REGION_2_LOOKUP}.arg

done

# Notes for installation of angsd v.0.921 using conda:
# ====================================================
# module load Anaconda2
# mkdir $DATA/angsd-env
# conda create --prefix $DATA/angsd-env --copy python=2.7
# source activate $DATA/angsd-env
# conda install -c bioconda angsd=0.921
# ====================================================
