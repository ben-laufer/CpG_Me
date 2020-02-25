#!/bin/bash
#
#SBATCH --job-name=CpG_Me_PE
#SBATCH --ntasks=1 # Number of cores/threads
#SBATCH --mem-per-cpu=1000 # Ram in Mb
#SBATCH --partition=production  
#SBATCH --output=CpG_Me_PE_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=CpG_Me_PE_%A_%a.err # File to which STDERR will be written
#SBATCH --time=0-00:10:00

##########################################################################################
# Author: Ben Laufer
# Email: blaufer@ucdavis.edu 
##########################################################################################

##############
# Initialize #
##############

# Manually set mainPath and partition

export mainPath="/share/lasallelab"
export partition="production"

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

jid1=$(sbatch \
--job-name=Trim \
--output=Trim_%j.out \
--error=Trim_%j.err \
--partition=${partition} \
--ntasks=15 \
--mem=12000 \
--time=0-03:00:00 \
${mainPath}/programs/CpG_Me/Paired-end/CpG_Me_PE_switch.sh \
trim \
${genome} \
| cut -d " " -f 4)

####################
# Screen and Align #
####################

# Set threads for fastqscreen in config file
# Each multicore needs 3 cores and 5 GB RAM per a core for directional libraries

jid2=$(sbatch \
--job-name=Align \
--output=Align_%j.out \
--error=Align_%j.err \
--partition=${partition} \
--dependency=afterok:${jid1} \
--ntasks=18 \
--mem=64000 \
--time=3-00:00:00 \
${mainPath}/programs/CpG_Me/Paired-end/CpG_Me_PE_switch.sh \
align \
${genome} \
| cut -d " " -f 4)

#########################
# Remove PCR Duplicates #
#########################

jid3=$(sbatch \
--job-name=Dedup \
--output=Dedup_%j.out \
--error=Dedup_%j.err \
--partition=${partition} \
--dependency=afterok:${jid2} \
--ntasks=1 \
--mem=30000 \
--time=2-00:00:00 \
${mainPath}/programs/CpG_Me/Paired-end/CpG_Me_PE_switch.sh \
deduplicate \
${genome} \
| cut -d " " -f 4) 

#######################
# Insert Size Metrics #
#######################

jid4=$(sbatch \
--job-name=Insert \
--output=Insert_%j.out \
--error=Insert_%j.err \
--partition=${partition} \
--dependency=afterok:${jid3} \
--ntasks=2 \
--mem=10000 \
--time=0-04:00:00 \
${mainPath}/programs/CpG_Me/Paired-end/CpG_Me_PE_switch.sh \
insert \
| cut -d " " -f 4) 

#######################
# Nucleotide Coverage #
#######################

jid5=$(sbatch \
--job-name=Coverage \
--output=Coverage_%j.out \
--error=Coverage_%j.err \
--partition=${partition} \
--dependency=afterok:${jid3} \
--ntasks=1 \
--mem=4000 \
--time=2-00:00:00 \
${mainPath}/programs/CpG_Me/Paired-end/CpG_Me_PE_switch.sh \
coverage \
${genome} \
| cut -d " " -f 4) 

#######################
# Extract Methylation #
#######################

# Each multicore needs 3 cores, 2GB overhead on buffer --split_by_chromosome \
jid6=$(sbatch \
--job-name=Extract \
--output=Extract_%j.out \
--error=Extract_%j.err \
--partition=${partition} \
--dependency=afterok:${jid3} \
--ntasks=18 \
--mem-per-cpu=2000 \
--time=2-00:00:00 \
${mainPath}/programs/CpG_Me/Paired-end/CpG_Me_PE_switch.sh \
extract \
${genome} \
| cut -d " " -f 4)

###########################
# Cytosine and QC reports #
###########################

sbatch \
--job-name=Report \
--output=Report_%j.out \
--error=Report_%j.err \
--partition=${partition} \
--dependency=afterok:${jid5}:${jid6} \
--ntasks=3 \
--mem-per-cpu=2000 \
--time=2-00:00:00 \
${mainPath}/programs/CpG_Me/Paired-end/CpG_Me_PE_switch.sh \
cytosineReport \
${genome}

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo ${runtime}

squeue -u $USER -o "%.8A %.4C %.10m %.20E"
