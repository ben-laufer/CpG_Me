#!/bin/bash
#
#SBATCH --job-name=rn6
#SBATCH --workdir /share/lasallelab/genomes/
#SBATCH --partition=production               
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=4000
#SBATCH --output=rn6.out # File to which STDOUT will be written
#SBATCH --error=rn6.err # File to which STDERR will be written
#SBATCH --time=1-00:00:00

export mainPath="/share/lasallelab"
PATH="$PATH:${mainPath}/programs/CpG_Me/Bismark-master/"
module load bowtie2/2.3.4.1
module load samtools/1.9

bowtie2-build rn6/rn6.fa rn6
bismark_genome_preparation --bowtie2 --verbose ${mainPath}/genomes/rn6/
