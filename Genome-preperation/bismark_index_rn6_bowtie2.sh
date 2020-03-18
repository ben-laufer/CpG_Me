#!/bin/bash
#
#SBATCH --job-name=rn6
#SBATCH --workdir /share/lasallelab/genomes/
#SBATCH --partition=production               
#SBATCH --ntasks=12
#SBATCH --mem=18000
#SBATCH --output=rn6.out # File to which STDOUT will be written
#SBATCH --error=rn6.err # File to which STDERR will be written
#SBATCH --time=0-05:00:00

export mainPath="/share/lasallelab"
PATH="$PATH:${mainPath}/programs/CpG_Me/Bismark-master/"
module load bowtie2/2.3.4.1
module load samtools/1.10

mkdir ${mainPath}/genomes/rn6
cd ${mainPath}/genomes/rn6

rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/895/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_genomic.fna.gz .

gunzip GCF_000001895.5_Rnor_6.0_genomic.fna.gz

mv GCF_000001895.5_Rnor_6.0_genomic.fna.gz rn6.fa

bowtie2-build --threads 12 rn6.fa rn6

mkdir Bowtie2
mv *.bt2 Bowtie2/

bismark_genome_preparation --bowtie2 --parallel 6 --verbose .
