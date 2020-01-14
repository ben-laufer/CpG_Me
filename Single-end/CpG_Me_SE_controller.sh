#!/bin/bash
#
#SBATCH --job-name=CpG_Me_SE
#SBATCH --ntasks=1 # Number of cores/threads
#SBATCH --mem-per-cpu=1000 # Ram in Mb
#SBATCH --partition=production 
#SBATCH --output=CpG_Me_SE_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=CpG_Me_SE_%A_%a.err # File to which STDERR will be written
#SBATCH --time=0-00:10:00

##########################################################################################
# Author: Ben Laufer
# Email: blaufer@ucdavis.edu 
##########################################################################################

##############
# Initialize #
##############

# Manually set mainPath

export mainPath="/share/lasallelab"

# Command line arguments set genome and array variables
# Provide a task_samples.txt file of sample ids (no file extensions) with one per a line in working directory and a raw_sequences folder with paired fastq files (.fq.gz)

genome=$1  

###################
# Run Information #
###################

start=`date +%s`

hostname

########
# Trim #
########

# M-bias correction
jid1=$(sbatch \
--ntasks=3 \
--mem-per-cpu=1000 \
--time=1-00:00:00 \
${mainPath}/programs/CpG_Me/Single-end/CpG_Me_SE_switch.sh \
trim \
${genome} \
| cut -d " " -f 4)

####################
# Screen and Align #
####################

# Set threads for fastqscreen in config file
# Each multicore needs 3 cores and 5 GB RAM per a core for directional libraries

jid2=$(sbatch \
--dependency=afterok:$jid1 \
--ntasks=18 \
--mem-per-cpu=5000 \
--time=5-00:00:00 \
${mainPath}/programs/CpG_Me/Single-end/CpG_Me_SE_switch.sh \
align \
${genome} \
| cut -d " " -f 4)

#########################
# Remove PCR Duplicates #
#########################

jid3=$(sbatch \
--dependency=afterok:$jid2 \
--ntasks=1 \
--mem-per-cpu=30000 \
--time=1-12:00:00 \
${mainPath}/programs/CpG_Me/Single-end/CpG_Me_SE_switch.sh \
deduplicate \
${genome} \
| cut -d " " -f 4) 

#######################
# Nucleotide Coverage #
#######################

jid4=$(sbatch \
--dependency=afterok:$jid3 \
--ntasks=1 \
--mem-per-cpu=4000 \
--time=2-00:00:00 \
${mainPath}/programs/CpG_Me/Single-end/CpG_Me_SE_switch.sh \
coverage \
${genome} \
| cut -d " " -f 4) 

#######################
# Extract Methylation #
#######################

# Each multicore needs 3 cores, 2GB overhead on buffer --split_by_chromosome \
jid5=$(sbatch \
--dependency=afterok:$jid3 \
--ntasks=18 \
--mem-per-cpu=2000 \
--time=2-00:00:00 \
${mainPath}/programs/CpG_Me/Single-end/CpG_Me_SE_switch.sh \
extract \
${genome} \
| cut -d " " -f 4)

############################
# Merge CpGs and QC report #
############################

# Generate merged CpG methylation for bsseq DMRfinder 
# Merge CpGs is an experimental feature
jid6=$(sbatch \
--dependency=afterok:$jid5 \
--ntasks=3 \
--mem-per-cpu=2000 \
--time=2-00:00:00 \
${mainPath}/programs/CpG_Me/Single-end/CpG_Me_SE_switch.sh \
mergeCpGs \
${genome} \
| cut -d " " -f 4)

###########################################
# DSS/DMRfinder and WGBS_tools Conversion #
###########################################

sbatch \
--dependency=afterok:$jid6 \
--ntasks=1 \
--mem-per-cpu=25000 \
--time=0-00:20:00 \
${mainPath}/programs/CpG_Me/Single-end/CpG_Me_SE_switch.sh \
format \
${genome}

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo $runtime

squeue -u $USER -o "%.8A %.4C %.10m %.20E"
