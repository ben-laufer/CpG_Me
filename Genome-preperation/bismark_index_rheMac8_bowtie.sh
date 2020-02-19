#!/bin/bash
#
#SBATCH --job-name=rheMac8
#SBATCH --workdir /share/lasallelab/genomes/
#SBATCH --partition=production               
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=4000
#SBATCH --output=rheMac8.out # File to which STDOUT will be written
#SBATCH --error=rheMac8.err # File to which STDERR will be written
#SBATCH --time=1-00:00:00

export mainPath="/share/lasallelab"
PATH="$PATH:${mainPath}/programs/CpG_Me/Bismark-master/"
module load bowtie2/2.3.4.1
module load samtools/1.9

bowtie2-build rheMac8/rheMac8.fa rheMac8
bismark_genome_preparation --bowtie2 --verbose ${mainPath}/genomes/rheMac8/
