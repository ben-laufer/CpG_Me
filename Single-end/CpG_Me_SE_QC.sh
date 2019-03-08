#!/bin/bash
#
#SBATCH --job-name=CpG_Me_SE
#SBATCH --ntasks=1 # Number of cores/threads
#SBATCH --mem-per-cpu=2000 # Ram in Mb
#SBATCH --partition=production    
#SBATCH --output=CpG_Me_SE_QC_%A.out # File to which STDOUT will be written
#SBATCH --error=CpG_Me_SE_QC_%A.err # File to which STDERR will be written
#SBATCH --time=0-06:00:00

##########################################################################################
# Author: Ben Laufer
# Email: blaufer@ucdavis.edu 
##########################################################################################

##############
# Initialize #
##############

# Manually set mainPath

mainPath=/share/lasallelab

###################
# Run Information #
###################

start=`date +%s`

hostname
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "My SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

THREADS=${SLURM_NTASKS}
MEM=$(expr ${SLURM_MEM_PER_CPU} / 1024)

echo "Allocated threads: " $THREADS
echo "Allocated memory: " $MEM

################
# Load Modules #
################

module load bismark/0.20.1
module load multiqc/1.7

###########
# MultiQC #
###########

call="multiqc
. \
 --config ${mainPath}/programs/CpG_Me/Single-end/multiqc_config_SE.yaml"

echo $call
eval $call

###########
# Bismark #
###########

call="bismark2summary \
"$(find `.` -name '*_bt2.bam' -print | tr '\n' ' ')""

echo $call
eval $call

#########
# Tidy  #
#########

# Remove non-deduplicated BAM files
if [ -f "bismark_summary_report.html" ] ; then
    find . -type f -name "*_bt2.bam" -exec rm -f {} \;
fi

# Copy cytosine reports to central directory
mkdir cytosine_reports
"$(find `.` -name '*cov.gz.CpG_report.txt.gz' -print0 | xargs -0 cp -t cytosine_reports)"

# Copy merged cytosine reports to central directory
mkdir cytosine_reports_merged
"$(find `.` -name '*merged_CpG_evidence.cov.gz' -print0 | xargs -0 cp -t cytosine_reports_merged)"
 
###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo $runtime
