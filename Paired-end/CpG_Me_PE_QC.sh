#!/bin/bash
#
#SBATCH --job-name=CpG_QC
#SBATCH --ntasks=2 # Number of cores/threads
#SBATCH --mem-per-cpu=4000 # Ram in Mb
#SBATCH --partition=production  
#SBATCH --output=CpG_Me_PE_QC_%j.out # File to which STDOUT will be written
#SBATCH --error=CpG_Me_PE_QC_%j.err # File to which STDERR will be written
#SBATCH --time=0-6:00:00

##########################################################################################
# Author: Ben Laufer
# Email: blaufer@ucdavis.edu 
##########################################################################################

##############
# Initialize #
##############

# Manually set mainPath

export mainPath="/share/lasallelab"

###################
# Run Information #
###################

start=`date +%s`

hostname

THREADS=${SLURM_NTASKS}
MEM=$((${SLURM_MEM_PER_CPU}/1024))

echo "Allocated threads: ${THREADS}"
echo "Allocated memory:  ${MEM}"

################
# Load Modules #
################

PATH="$PATH:${mainPath}/programs/CpG_Me/Bismark-master/"
module load bowtie2/2.3.4.1
module load samtools/1.10
module load multiqc/1.8

#########
# Tidy  #
#########

mkdir slurm_logs
mv {*.out,*.err} ./slurm_logs

###########
# MultiQC #
###########

call="multiqc \
. \
--ignore slurm_logs/ \
--ignore raw_sequences/ \
--config ${mainPath}/programs/CpG_Me/Paired-end/multiqc_config_PE.yaml"

echo ${call}
eval ${call}

###########
# Bismark #
###########

call="bismark2summary \
"$(find `.` -name '*_pe.bam' -print | tr '\n' ' ')""

echo ${call}
eval ${call}

#########
# Tidy  #
#########

# Remove non-deduplicated BAM files
if [ -f "bismark_summary_report.html" ]
then
    find . -type f -name "*_pe.bam" -exec rm -f {} \;
fi

# Copy cytosine reports to central directory
mkdir cytosine_reports
find \
. \
-name '*cov.gz.CpG_report.txt.gz' \
-type f \
-not -path "./cytosine_reports" \
-print0 | \
xargs -0 \
cp -t \
"./cytosine_reports" 

# Copy merged cytosine reports to central directory
mkdir cytosine_reports_merged
find \
. \
-name '*merged_CpG_evidence.cov.gz' \
-type f \
-not -path "./cytosine_reports_merged" \
-print0 | \
xargs -0 \
cp -t \
"./cytosine_reports_merged" 

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo ${runtime}
