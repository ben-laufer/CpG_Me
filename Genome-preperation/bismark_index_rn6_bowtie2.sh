#!/bin/bash
#
#SBATCH --job-name=rn6
#SBATCH --workdir /share/lasallelab/genomes/
#SBATCH --partition=production               
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=4000
#SBATCH --output=slurmlogs/rn6_%A.out # File to which STDOUT will be written
#SBATCH --error=slurmlogs/rn6_%A.err # File to which STDERR will be written
#SBATCH --time=1-00:00:00

module load bowtie2/2.3.4.1
module load samtools/1.8
module load bismark/0.20.0

bowtie2-build rn6/rn6.fa rn6
bismark_genome_preparation --bowtie2 --verbose /share/lasallelab/genomes/rn6/
