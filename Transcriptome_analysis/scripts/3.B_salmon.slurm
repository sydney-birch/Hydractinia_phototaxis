#!/bin/bash

## Note - Slurm script comments require two hash symbols (##).  A single
## hash symbol immediately followed by SBATCH indicates an SBATCH
## directive.  "##SBATCH" indicates the SBATCH command is commented
## out and is inactive.

## NTasks is not thread count, be sure to leave it set at 1
#SBATCH --ntasks=1

## If your program will be using less than 24 threads, or you
## require more than 24 threads, set cpus-per-task to the 
## desired threadcount.  Leave this commented out for the
## default 24 threads.
#SBATCH --cpus-per-task=24

## You will need to specify a minimum amount of memory in the
## following situaitons:
##   1. If you require more than 128GB of RAM, specify either:
##      a. "--mem=512000" for at least 512GB of RAM (6 possible nodes)
##      b. "--mem=1000000" for at least 1TB of RAM (2 possible nodes)
##   2. If you are running a job with less than 24 threads, you will
##      normally be given your thread count times 5.3GB in RAM.  So
##      a single thread would be given about 5GB of RAM.  If you
##      require more, please specify it as a "--mem=XXXX" option,
##      but avoid using all available RAM so others may share the node.
##SBATCH --mem=512000

## Normally jobs will restart automatically if the cluster experiences
## an unforeseen issue.  This may not be desired if you want to retain
## the work that's been performed by your script so far.   
## --no-requeue

## Normal Slurm options
## SBATCH -p shared
#SBATCH --job-name="salmon"
#SBATCH --output=salmon_orp.output

## Load the appropriate modules first.  Linuxbrew/colsa contains most
## programs, though some are contained within the anaconda/colsa
## module.  Refer to http://premise.sr.unh.edu for more info.
module purge
module load linuxbrew/colsa


## Instruct your program to make use of the number of desired threads.
## As your job will be allocated an entire node, this should normally
## be 24.

##format actinulae ORP reference transcriptome - make index 	
salmon index -t actinula_total.ORP.fa-mod.fa -i actinulae_orp_index --type quasi -k 31  



##map libraries 
##STG_1_R1
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_1_R1.R1.fastq.gz -2 raw_larva_data/STG_1_R1.R2.fastq.gz -o STG_1_R1 

##STG_1_R2
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_1_R2.R1.fastq.gz -2 raw_larva_data/STG_1_R2.R2.fastq.gz -o STG_1_R2 

##STG_1_R3
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_1_R3.R1.fastq.gz -2 raw_larva_data/STG_1_R3.R2.fastq.gz -o STG_1_R3 

##STG_1_R4
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_1_R4.R1.fastq.gz -2 raw_larva_data/STG_1_R4.R2.fastq.gz -o STG_1_R4 

##STG_1_R5
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_1_R5.R1.fastq.gz -2 raw_larva_data/STG_1_R5.R2.fastq.gz -o STG_1_R5 

##STG_1_R6
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_1_R6.R1.fastq.gz -2 raw_larva_data/STG_1_R6.R2.fastq.gz -o STG_1_R6 


##STG_2_R1
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_2_R1.R1.fastq.gz -2 raw_larva_data/STG_2_R1.R2.fastq.gz -o STG_2_R1 

##STG_2_R2
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_2_R2.R1.fastq.gz -2 raw_larva_data/STG_2_R2.R2.fastq.gz -o STG_2_R2 

##STG_2_R3
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_2_R3.R1.fastq.gz -2 raw_larva_data/STG_2_R3.R2.fastq.gz -o STG_2_R3 

##STG_2_R4
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_2_R4.R1.fastq.gz -2 raw_larva_data/STG_2_R4.R2.fastq.gz -o STG_2_R4 

##STG_2_R5
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_2_R5.R1.fastq.gz -2 raw_larva_data/STG_2_R5.R2.fastq.gz -o STG_2_R5 

##STG_2_R6
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_2_R6.R1.fastq.gz -2 raw_larva_data/STG_2_R6.R2.fastq.gz -o STG_2_R6 


##STG_3_R1
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_3_R1.R1.fastq.gz -2 raw_larva_data/STG_3_R1.R2.fastq.gz -o STG_3_R1 

##STG_3_R2
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_3_R2.R1.fastq.gz -2 raw_larva_data/STG_3_R2.R2.fastq.gz -o STG_3_R2 

##STG_3_R3
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_3_R3.R1.fastq.gz -2 raw_larva_data/STG_3_R3.R2.fastq.gz -o STG_3_R3 

##STG_3_R4
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_3_R4.R1.fastq.gz -2 raw_larva_data/STG_3_R4.R2.fastq.gz -o STG_3_R4 

##STG_3_R5
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_3_R5.R1.fastq.gz -2 raw_larva_data/STG_3_R5.R2.fastq.gz -o STG_3_R5 

##STG_3_R6
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_3_R6.R1.fastq.gz -2 raw_larva_data/STG_3_R6.R2.fastq.gz -o STG_3_R6 



##STG_4_R1
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_4_R1.R1.fastq.gz -2 raw_larva_data/STG_4_R1.R2.fastq.gz -o STG_4_R1 

##STG_4_R2
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_4_R2.R1.fastq.gz -2 raw_larva_data/STG_4_R2.R2.fastq.gz -o STG_4_R2 

##STG_4_R3
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_4_R3.R1.fastq.gz -2 raw_larva_data/STG_4_R3.R2.fastq.gz -o STG_4_R3 

##STG_4_R4
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_4_R4.R1.fastq.gz -2 raw_larva_data/STG_4_R4.R2.fastq.gz -o STG_4_R4 

##STG_4_R5
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_4_R5.R1.fastq.gz -2 raw_larva_data/STG_4_R5.R2.fastq.gz -o STG_4_R5 

##STG_4_R6
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_4_R6.R1.fastq.gz -2 raw_larva_data/STG_4_R6.R2.fastq.gz -o STG_4_R6 



##STG_5_R1
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_5_R1.R1.fastq.gz -2 raw_larva_data/STG_5_R1.R2.fastq.gz -o STG_5_R1 

##STG_5_R2
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_5_R2.R1.fastq.gz -2 raw_larva_data/STG_5_R2.R2.fastq.gz -o STG_5_R2 

##STG_5_R3
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_5_R3.R1.fastq.gz -2 raw_larva_data/STG_5_R3.R2.fastq.gz -o STG_5_R3 

##STG_5_R4
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_5_R4.R1.fastq.gz -2 raw_larva_data/STG_5_R4.R2.fastq.gz -o STG_5_R4 

##STG_5_R5
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_5_R5.R1.fastq.gz -2 raw_larva_data/STG_5_R5.R2.fastq.gz -o STG_5_R5 

##STG_5_R6
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_5_R6.R1.fastq.gz -2 raw_larva_data/STG_5_R6.R2.fastq.gz -o STG_5_R6 




##STG_6_R1
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_6_R1.R1.fastq.gz -2 raw_larva_data/STG_6_R1.R2.fastq.gz -o STG_6_R1 

##STG_6_R2
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_6_R2.R1.fastq.gz -2 raw_larva_data/STG_6_R2.R2.fastq.gz -o STG_6_R2 

##STG_6_R3
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_6_R3.R1.fastq.gz -2 raw_larva_data/STG_6_R3.R2.fastq.gz -o STG_6_R3 

##STG_6_R4
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_6_R4.R1.fastq.gz -2 raw_larva_data/STG_6_R4.R2.fastq.gz -o STG_6_R4 

##STG_6_R5
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_6_R5.R1.fastq.gz -2 raw_larva_data/STG_6_R5.R2.fastq.gz -o STG_6_R5 

##STG_6_R6
salmon quant -p 8 --seqBias --gcBias -i actinulae_orp_index -l a -1 raw_larva_data/STG_6_R6.R1.fastq.gz -2 raw_larva_data/STG_6_R6.R2.fastq.gz -o STG_6_R6 

