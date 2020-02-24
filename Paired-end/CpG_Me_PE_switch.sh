#!/bin/bash

##########################################################################################
# Author: Ben Laufer
# Email: blaufer@ucdavis.edu 
##########################################################################################

##############
# Initialize #
##############

# Manually set mainPath

export mainPath="/share/lasallelab"

# Command line arguments set the module and genome variables
# Provide a task_samples.txt file of sample ids (no file extensions) with one per a line in working directory and a raw_sequences folder with paired fastq files (.fq.gz)

module=$1  
genome=$2     

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

module load trim_galore/0.6.5
PATH="$PATH:${mainPath}/programs/CpG_Me/Bismark-master/"
module load bowtie2/2.3.4.1
module load samtools/1.9
PATH="$PATH:${mainPath}/programs/CpG_Me/fastq_screen_v0.14.0/"
export PYTHON_EGG_CACHE="${mainPath}/programs/CpG_Me"
module load picard-tools/2.18.4

######################
# Set Up Environment #
######################

directory=${PWD}/
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" task_samples.txt`
rawpath=${directory}raw_sequences/
mappath=${directory}${sample}
fastq1=${rawpath}${sample}_1.fq.gz
fastq2=${rawpath}${sample}_2.fq.gz
output=${mappath}${sample}
trim1=${sample}_1_val_1.fq.gz
trim2=${sample}_2_val_2.fq.gz
BAM=${sample}_1_val_1_bismark_bt2_pe.bam
dedupBAM=${sample}_1_val_1_bismark_bt2_pe.deduplicated.bam
sortedBAM=${sample}_1_val_1_bismark_bt2_pe.deduplicated.sorted.bam
insert=${sample}_1_val_1_bismark_bt2_pe.deduplicated.sorted.bam.insert.txt
histogram=${sample}_1_val_1_bismark_bt2_pe.deduplicated.sorted.bam.histogram.pdf
cov=${sample}_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
CpH=Non_CpG_context_${sample}_1_val_1_bismark_bt2_pe.deduplicated.txt.gz
CpGmerge=${sample}_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.merged_CpG_evidence.cov.gz

#################
# Case Switches #
#################

case ${module} in
     trim)      
          ########
          # Trim #
          ########
          
          mkdir ${mappath}

 		  # Use 2color for NovaSeq and NextSeq, replace with quality for HiSeq and MiSeq
          # M-bias correction Swift's Accel NGS 
          call="trim_galore \
          --paired \
          --cores 4 \
          --2colour 20 \
          --fastqc \
          --clip_r1 10 \
          --clip_r2 20 \
          --three_prime_clip_r1 10 \
          --three_prime_clip_r2 10 \
          --output_dir ${mappath} \
          ${fastq1} \
          ${fastq2}" 

          echo ${call}
          eval ${call}
          ;;     
     align)      
          ##########
          # Screen #
          ##########
          
          cd ${mappath}

          call="fastq_screen \
          --conf ${mainPath}/programs/CpG_Me/fastq_screen.conf \
          --bisulfite \
          ${trim1} \
          ${trim2}" 

          echo ${call}
          eval ${call}

          #########
          # Align #
          #########

          cd ${mappath}

          # Each multicore needs 3 cores and 5 GB RAM per a core for directional libraries

          call="bismark \
          -n 1 \
          --genome ${mainPath}/genomes/${genome}/ \
          --multicore 6 \
          -1 ${trim1} \
          -2 ${trim2}"

          echo ${call}
          eval ${call}
              
          #############################
          # Remove Intermediate Files #
          #############################

          if [ -f ${BAM} ]
          then
            rm ${trim1}
          	rm ${trim2}
          fi
          ;;     
     deduplicate)
          #########################
          # Remove PCR Duplicates #
          #########################

          cd ${mappath}

          # WGBS only
          call="deduplicate_bismark \
          --bam \
          --paired \
          ${BAM}"

          echo ${call}
          eval ${call}
          ;;
     insert)
     	  #######################
     	  # Insert Size Metrics #
     	  #######################
     	  
     	  cd ${mappath}
     	  
     	  # Can't pipe because of multiQC thinking all are named stdin
     	  call="picard SortSam \
		  INPUT=${dedupBAM} \
		  OUTPUT=${sortedBAM} \
		  SORT_ORDER=coordinate"
		  
		  echo ${call}
		  eval ${call}
		  
		  call="picard CollectInsertSizeMetrics \
		  INPUT=${sortedBAM} \
		  OUTPUT=${insert} \
		  HISTOGRAM_FILE=${histogram} \
		  ASSUME_SORTED=TRUE"

		  echo ${call}
		  eval ${call}
		  
		  rm ${sortedBAM}
     	  ;;
     coverage)
          #######################
          # Nucleotide Coverage #
          #######################
          
          cd ${mappath}

          call="bam2nuc \
          --genome_folder ${mainPath}/genomes/${genome}/ \
          ${dedupBAM}"

          echo ${call}
          eval ${call}
          ;; 
     extract)
          #######################
          # Extract Methylation #
          #######################

          cd ${mappath}

          # Each multicore needs 3 cores, 2GB overhead on buffer --split_by_chromosome \
          # Use --scaffolds for genomes with many contigs
          if [ ${genome} == "rheMac8" ]
          then
          	call="bismark_methylation_extractor \
          	--paired-end \
          	--gzip \
          	--comprehensive \
          	--merge_non_CpG \
          	--bedGraph \
          	--multicore 6 \
          	--scaffolds \
          	--buffer_size 34G \
          	${dedupBAM}"

          	echo ${call}
          	eval ${call}
          
          else
          	call="bismark_methylation_extractor \
          	--paired-end \
          	--gzip \
          	--comprehensive \
          	--merge_non_CpG \
          	--bedGraph \
          	--multicore 6 \
          	--buffer_size 34G \
          	${dedupBAM}"

          	echo ${call}
          	eval ${call}	
          fi
          
          rm ${CpH}
          ;;
     cytosineReport)
          ###################
          # Cytosine Report #
          ###################

          cd ${mappath}

          # Generate genome-wide cytosine reports
          # Merge CpGs is an experimental feature, so DMRichR doesn't use those
          call="coverage2cytosine \
          --output ${mappath}/${cov} \
          --genome_folder ${mainPath}/genomes/${genome}/ \
          --gzip \
          --merge_CpG \
          ${mappath}/${cov}"

          echo ${call}
          eval ${call}

          #############
          # QC Report #
          #############

          call="bismark2report"

          echo ${call}
          eval ${call}
          ;;
     *)
          echo "Error: Pipeline case selection of ${module} is invalid or not specified."
          echo "Please select either trim, align, deduplicate, insert, coverage, extract, or cytosineReport"
          ;;
esac

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo ${runtime}
