#!/bin/bash
#
#SBATCH --job-name=mm10
#SBATCH --workdir /share/lasallelab/genomes/
#SBATCH --partition=production               
#SBATCH --ntasks=12
#SBATCH --mem=18000
#SBATCH --output=mm10.out # File to which STDOUT will be written
#SBATCH --error=mm10.err # File to which STDERR will be written
#SBATCH --time=0-05:00:00

export mainPath="/share/lasallelab"
PATH="$PATH:${mainPath}/programs/CpG_Me/Bismark-master/"
module load bowtie2/2.3.4.1
module load samtools/1.10

mkdir ${mainPath}/genomes/mm10
cd ${mainPath}/genomes/mm10

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz .
gunzip mm10.fa.gz

bowtie2-build --threads 12 --verbose mm10.fa mm10

mkdir Bowtie2
mv *.bt2 Bowtie2/

bismark_genome_preparation --bowtie2 --parallel 6 --verbose .
