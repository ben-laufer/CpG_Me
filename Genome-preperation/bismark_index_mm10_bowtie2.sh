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

rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz .

gunzip GCF_000001635.26_GRCm38.p6_genomic.fna.gz

mv GCF_000001635.26_GRCm38.p6_genomic.fna mm10.fa

bowtie2-build --threads 12 mm10.fa mm10

mkdir Bowtie2
mv *.bt2 Bowtie2/

bismark_genome_preparation --bowtie2 --parallel 6 --verbose .
