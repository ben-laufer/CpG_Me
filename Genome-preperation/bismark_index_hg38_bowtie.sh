#!/bin/bash
#
#SBATCH --job-name=hg38
#SBATCH --workdir /share/lasallelab/genomes/
#SBATCH --partition=production               
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=4000
#SBATCH --output=hg38.out # File to which STDOUT will be written
#SBATCH --error=hg38.err # File to which STDERR will be written
#SBATCH --time=1-00:00:00

export mainPath="/share/lasallelab"
PATH="$PATH:${mainPath}/programs/CpG_Me/Bismark-master/"
module load bowtie2/2.3.4.1
module load samtools/1.9

bowtie2-build hg38/hg38.fa hg38
bismark_genome_preparation --bowtie2 --verbose ${mainPath}/genomes/hg38/
