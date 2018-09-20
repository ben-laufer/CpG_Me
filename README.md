# CpG_Me
### A whole-genome bisulfite sequencing (WGBS) pipeline for the analysis of DNA methylation

CpG_Me_PE is a series of shell scripts that automate a WGBS workflow that takes you from raw fastq files to extracted CpG methylation count values, where it preprocesses data to remove biases and provides ample QC/QA. 

## Chastity Filtering

This workflow assumes your data is Illumina quality/chastity filtered, which most service providers these days will do by default.

You can check by using the following command, where file.fastq.gz represents your file:

`zcat JLBL001.fastq.gz | head -n 50`

If they aren’t you can accomplish this on command line via, where you change JLBL001 to your sample name

`zcat JLBL001*fastq.gz | zgrep -A 3 '^@.* [^:]*:N:[^:]*:' | zgrep -v "^--$" | gzip > JLBL001_filtered.fq.gz`

## How to use CpG_Me:
1.	Create a parent directory 
2.	Within that parent directory, add a text file called “task_samples.txt”, where each new line contains the entire sample name exactly as it appears on the fastq read pair files, aside from the end part (“_1.fq.gz” or “_2.fq.gz”). Only name a sample once, NOT twice, and make sure it is .fq.gz and not fastq.gz. Also, if you’re using excel or a windows desktop, you will need to change the linebreaks from windows to unix, which can be done using text wrangler.
3.	Within that parent directory create a folder called “raw_sequences” that contains all raw paired fastq files (.fq.gz)
4.	Now it’s ready to run, so FROM the parent directory, modify and run this command:

`sbatch --array=1-12 /share/lasallelab/programs/CpG_Me/CpG_Me_PE_controller.sh  hg38`

Let’s break this apart:
1)	sbatch is how you submit a job to a HPCC with a slurm workload manager
2)	--array=12 lets you specify the number of samples, as well as subset. Here we are running samples 1 to 12. You could run select samples using the following format --array=2,4-12
3)	The next call is the location of the executable shell script that will schedule all jobs with proper resources and dependencies on a per sample basis
4)	Genome (hg38, rheMac8, mm10, rn6)

There is also a final QC report to be run AFTER all samples have finished, which you also need to launch from the working directory

`sbatch /share/lasallelab/programs/CpG_Me/CpG_Me_QC_PE.sh` 
