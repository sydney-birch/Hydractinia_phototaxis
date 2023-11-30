# Step by Step Methods of Hydractinia Planula Transcriptome Analysis and Outputs

This file is the step by step instructions of our Transcriptome analysis. You will find the scripts used and an explanation of what they do as well as outputs. I am following the Transcriptome workflow we previously developed, which can be found in my Actinula_paper repository.

### 1. Prep Reads for Transcriptome Assembly
We sequenced planulae larvae over 4 developmental days: Day 1, Day 2, Day 3 (when known to be phototactic), and Day 4. We have 3 replicates for each of the 4 developmental Days (we will call stages here). The first step of our developmental transcriptome assembly is to prepare our reads for assembling a reference transcriptome. Here in this step you can either run the full cohesive script: 1_Full_Ref_Transcriptome_prep.py or run it as 2 scripts (1_steps_1-3a_v3.py and 1_steps_3b-4b.py) but with running it as 2 scripts you have to hardcode a collection in the second script (I would recommend running just the cohesive script which is what I demonstrate below).

**The goal of 1_Full_Ref_Transcriptome_prep.py** is to identify the highest quality replicate for each of the 6 stages and then concatenate the high-quality reps into representitive R1 and R2 files to be used in the assembly. The way this script identifies high-quality reps is by running FastQC on each R1 and R2 file for all reps in all stages and then calculating the following basic stats from the FastQC data:

1. ID the number of nucleotides: num_nucs = Number of Sequences * Read Length.
2. ID Normalized Quality: norm_qual = (Quality (the avg of means from FastQC file) - total min across all files)/(total max - total min).
3. Calculate Final Score: final_score = num_nucs * norm_qual.

The script then compares each of the final scores within each stage and then concatenates the 6 highest quality reads (1 from each stage) into a total_R1.fastq.gz and total_R2.fastq.gz (the output files). The script will also report the time it took to run the full script and a FastQC dir with all the FastQC results (note: the FastQC step will take the longest).

To run this script you will need to organize your directory as such:
`mkdir ORP_Prep`.
`mkdir ORP_Prep/raw_larva_reads`.

Within raw_larva_reads dir you should have subdirectories for each sample. Rename your subdirs with this format: group_info_additional-info, for example: STG_1_R1, STG_1_R2, STG_1_R3, etc. If you need to add additional info after the R# add a dash and the info (STG_1_R3-additional-info). The group for this example will be STG1 and within that group there will be STG_1_R1, STG_1_R2, STG_1_R3, etc. This is important because the script will find the best representative sample within each of your groups which will be concatenated in the total_R1.fastq.gz and total_R2.fastq.gz files. For the actinula data there are a total of 36 sub dirs with a total of 72 files (an R1 and R2 read file for each sample). Side Note: You do not need to rename the actual sample files (the R1 and R1 files within each sub directory) and keep everything zipped.

##### Once your directory is set up, run the script:
`./1_Full_Ref_Transcriptome_prep.py -d raw_larva_reads`

The highest quality replicates that are used in the Reference transcriptome are:
winners_R1 = ['STG_1_DP-18.R1', 'STG_2_DP-20.R1', 'STG_4_DP-28.R1', 'STG_3_DP-26.R1']
winners_R2 = ['STG_1_DP-18.R2', 'STG_2_DP-20.R2', 'STG_4_DP-28.R2', 'STG_3_DP-26.R2']

### 2. Run Transcriptome Assembler: Oyster River Protocol (ORP)    
   Now that we have our representative R1 and R2 files, we can assemble the transcriptome. Here we used the ORP which generates 3 assemblies and merges them into 1 high quality assembly. More information on this method can be found here: https://oyster-river-protocol.readthedocs.io/en/latest/  
   
   ##### Submit the slurm:   
   `sbatch 2_ORP.slurm`  
   
   The code in this slurm (full slurm script can be found in scripts_for_analysis folder):   
   ```oyster.mk \   
   MEM=150 \   
   CPU=24 \   
   READ1=total_R1.fastq.gz \   
   READ2=total_R2.fastq.gz \   
   RUNOUT=hydractinia_total   
   ```

   Output Quality Metrics for the assembly:  
   
    QUALITY REPORT FOR: hydractinia_total using the ORP version 2.2.8
    THE ASSEMBLY CAN BE FOUND HERE: /mnt/oldhome/plachetzki/sjb1061/hydractinia_transcriptomics/ORP/assemblies/hydractinia_total.ORP.fasta 

    BUSCO SCORE ~~~~~~~~~~~~~~~~~~~~~~>      C:100.0%[S:48.2%,D:51.8%],F:0.0%,M:0.0%,n:303
    TRANSRATE SCORE ~~~~~~~~~~~~~~~~~~>      0.39257
    TRANSRATE OPTIMAL SCORE ~~~~~~~~~~>      0.444
    UNIQUE GENES ORP ~~~~~~~~~~~~~~~~~>      13186
    UNIQUE GENES TRINITY ~~~~~~~~~~~~~>      11449
    UNIQUE GENES SPADES55 ~~~~~~~~~~~~>      12060
    UNIQUE GENES SPADES75 ~~~~~~~~~~~~>      10900
    UNIQUE GENES TRANSABYSS ~~~~~~~~~~>      11098
    READS MAPPED AS PROPER PAIRS ~~~~~>      93.87% 	


### 3. Quantify Reads - Run Salmon
   Before we quantify our reads, we first need to change the headers in our new assembly. Some of the headers are very long and will be difficult to work with when we run Transdecoder in the next step so it is better to re-format the headers now. 
   
   ##### To reformat the headers, Run 3.A_rename_fa_headers.py.   
   `./3.A_rename_fa_headers.py -a hydractinia_total.ORP.fasta -b Hs_planula`   
   
   The second argument here is the base of the new header, this script counts each header in order so the first header will look like: Hs_planula_t.1   
  
   The output from this script will be the modified fasta file:  hydractinia_total.ORP.fasta-mod.fa
     
   Now that the headers have been reformated, we can use this assembly in Salmon which quantifies reads. Info on salmon can be found here: https://salmon.readthedocs.io/en/latest/salmon.html#using-salmon. The two inputs you will need for this program are the assembly and all of the raw reads. This slurm script (3.B_salmon.slurm) will map each replicate to the assembly (note: you only have to make the index once). 
   
   ##### Run the salmon Slurm.  
   `sbatch 3.B_salmon.slurm`  
   
   Output: You will get directories for each sample - here I secure copied all of my output directories to my desktop and placed them in a folder named mapping. This will be used with EdgeR to find DEGs and when making heatmaps in R towards the end of this workflow. 
  
### 4. Run TransDecoder 
   Next we need to translate our assembly into protien space so we can use our assembly in Orthofinder in the next step. We will also have to alter the headers again keeping only valuable information and then we will reduce the assembly by using cd-hit to get rid of duplicates. For more info on transdecoder go here: https://github.com/TransDecoder/TransDecoder/wiki.  
   
   ##### First, run Transdecoder:    
   `sbatch 4.A_transdecoder.slurm`   
   
   The code in the slurm is:  
   `TransDecoder.LongOrfs -t hydractinia_total.ORP.fasta-mod.fa`.    
   Use the longest_orfs.pep file in the transdecoder dir - copy the file and rename it: hydractinia_total_ORP_prot.fa   
   
   ##### Next, rename the headers in the prot fasta:   
   `./4.B_rename_prot_headers.py -a hydractinia_total_ORP_prot.fa`      
   
   The headers after running Transdecoder look like:   
    >Gene.1::Hs_planula_t.2::g.1::m.1 type:complete len:719 gc:universal Hs_planula_t.2:105-2261(+)  
    	
   This script will take the last segment of the transdecoder header (>Gene.1::Hs_planula_t.2::g.1::m.1 type:complete len:719 gc:universal **Hs_planula_t.2:105-2261(+)**) and will remove the direction ((-) or (+)) and the colon to make the new header which would be: >Hs_planula_t.2..105-2261. This header provides the original transcript from which it came (t.#) and also provides where in that sequence it came (#-#). The output file will have the -mod.fa tag at the end of the input fasta file.   
   
   The output file is: hydractinia_total_ORP_prot.fa-mod.fa.  
   
   ##### Now to get rid of potential duplicates in the assembly, run cd-hit on our protien assembly:   
   `sbatch 4.C_cdhit.slurm`
   
   The code in the slurm:  
   `cd-hit -i hydractinia_total_ORP_prot.fa-mod.fa -o hydractinia_total_ORP_prot.fa-mod_reduced.fa -c 0.98`  
   
   output file: hydractinia_total_ORP_prot.fa-mod_reduced.fa  
   
   The number of transcripts: `grep ">" hydractinia_total.ORP.fasta-mod.fa | wc -l`  = **92784 transcripts**.  
   The number of Protien models: `grep ">" hydractinia_total_ORP_prot.fa-mod.fa | wc -l` = **73297 prot mods**.  
   
   The number of Protien models after cd-hit: `grep ">" hydractinia_total_ORP_prot-mod_reduced.fa | wc -l` = **48629 prot mods**.  
   
  ##### Side Note: At this point there should now be 5 fasta files:   
 * The original ORP nuc fasta: *hydractinia_total.ORP.fa*.  
 * The modified header ORP nuc fasta: *hydractinia_total.ORP.fasta-mod.fa*.  
 * The original Transdecoder prot fasta: *hydractinia_total_ORP_prot.fa*    
 * The modified header prot fasta: *hydractinia_total_ORP_prot.fa-mod.fa*.  
 * The reduced (cd-hit) prot fasta (w/ mod headers): *hydractinia_total_ORP_prot-mod_reduced.fa*.  
   
   From here on, we will be using the reduced prot fasta: *hydractinia_total_ORP_prot-mod_reduced.fa*.  
   
### 5. Run OrthoFinder
   Now that we have our reduced Hydractinia planula prot models, we can use it with OrthoFinder to identify orthogroups among taxa, for more info on OrthoFinder go here: https://github.com/davidemms/OrthoFinder. For this run, we will use a total of 7 taxa:  
  * actinula (*E. crocea* larva).  
  * Hydractinia planula (this transcriptome).   
  * Hydractinia adult polyp.  
  * *Nematostella*  
  * *Hydra*  
  * *Drosophila*.  
  * *Homo sapiens*.  
  
  ##### Run OrthoFinder slurm:  
  `sbatch 5_orthofinder.slurm`  
  
  The code in the slurm:   
  `orthofinder.py -f ./fastas/ -t 24 -M msa -S diamond`   
  
   The ouput files you need to transfer to your desktop to be used in R in later steps are: 
  * Orthogroups.tsv (renamed it to: Orthogroups_5-11-20.tsv).  
  * Orthogroups.GeneCount.tsv (renamed it to: Orthogroups.GeneCount_5-11-20.tsv).  


### 6. Download Gene Sets  

   The goal of step 6 is to download the protien fastas for genes in different sensory gene sets currated by the Broad Institute. This section is split up into two major steps, A and B. Within the B step there are 4 sub steps that need to be performed. This method uses the python module bio entrez to download each sequence within a gene set and it double checks that the correct number of sequences have been downloaded.     
   
   ##### A. Select gene sets from the Broad Institute Gene Set Enrichment Analysis (GSEA) https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp    
   We used 3 gene sets for this paper:   
   [GO_SENSORY_PERCEPTION_OF_LIGHT_STIMULUS](https://www.gsea-msigdb.org/gsea/msigdb/cards/GO_SENSORY_PERCEPTION_OF_LIGHT_STIMULUS.html).  
   [GO_SENSORY_SYSTEM_DEVELOPMENT](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_SENSORY_SYSTEM_DEVELOPMENT.html).  
   [GO_TRANSCRIPTION_FACTOR_ACTIVITY](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOMF_DNA_BINDING_TRANSCRIPTION_FACTOR_BINDING.html).  

   For each gene set, click on the show members link and then copy all info into seperate excel files and save as csv files. Make sure you have the entrez id, gene symbol, and the description.   
   
   ##### B. Download Gene set sequences *(run on each gene set)* 
   
   ###### B.1 Import CSV files to terminal 
   You can either secure copy csv files and send to terminal or open a nano window and copy and paste. 
   
   ###### B.2 Run clean up script: 6.B.2_clean_up_csv.py
   This script will go through the csv file and will clean up the descriptions, it will remove the trailing .. and any brackets with text.   
   `./6.B.2_clean_up_csv.py -i sensory_percep_light_stim_geneset.txt`   
   
   output: sensory_percep_chem_light_geneset.txt-mod   
   
   ###### B.3 Download Sequences using NCBI Entrez database.  
   There are two options here, both do the same thing but because of the number of sequences in a gene set, there needs to be a pause in the connection to the database. **Choose the correct option based on the number of genes in the gene set.**   
   
   There are 3 parts to downloading sequences in these scripts. The first part performs the intial search of genes in *Homo sapiens* to get accession IDs. This will generate two ouput files:   
   summary_info.txt which contains the info for each gene symbol (we need the gene ID from this file).   
   gene_tables.txt which has all of the gene information (we need the protien accession IDs from this file for each symbol).   
   
   The second part creates a dictionary where the gene symbols are the keys and accession ids are the values - this is populated by going through the gene_tables.txt file and adding accession ids that have been verified and start with NP_ . If you are using option 1 because your gene set has 1-150 genes the script will move right into part 3, if you are using option 2 then the script will write a temporary file of this dicitonary to be used in part 3.  
   
   The third part of the script will use the dictionary just created to search the accession id that starts with NP for each gene symbol in the human database and will write out the protien sequence to a FASTA file. This part will also generate a file called gene_symbol_accid which will be tab delimited with 2 columns with gene symbols in the first and accession ids in the second. These are the two files you will need in the following steps, I recommend copying those 2 files into a new directory for step 7.   
   
   Option 1: **1-150** genes in gene set use 6.B.3a_get_entrez_fastas-v5.py   
   `./6.B.3a_get_entrez_fastas-v5.py -i sensory_percep_light_stim_geneset.txt-mod -o sensory_percep_chem_stim.fa`   

   Option 2: **150-500** genes in gene set use 2 scripts:  
   `./6.B.3b_1_split_get_entrez_fasta-v5.py -i sensory_percep_light_stim_geneset.txt-mod`     
   `./6.B.3b_2_split_get_entrez_fasta-v5.py -o sensory_percep_light_stim.fa`    

###### B.4 Check for missing seqs: 6.B.4_check_missing_seqs-v2.py   
   Make sure you have downloaded all of the genes in the gene set. This script will compare one of the output files (gene_symbol_accid) with the cleaned up csv file from step B.2. If there are missing/incorrect genes you will have to change them by hand by searching the entrez gene id in the NCBI gene database. Make sure to adjust both the FASTA file and the gene_symbol_accid file accordingly. At the end, copy your FASTA file and your gene_symbol_accid file into a new directory to be used in the next step.  
  
  Run script  
  `./6.B.4_check_missing_seqs-v2.py -a sensory_percep_light_stim_geneset.txt-mod -b gene_symbol_accid`   


### 7. Find Human Representative Sequences for Gene Sets in our Human protein FASTA  
  Since we have modified the headers of our Human protien FASTA, we want to identify the NCBI sequences from the gene sets in our Human FASTA. So we are going to BLAST the NCBI gene set sequences (query) to our Human FASTA with altered headers. We will only allow 1 hit per query sequence. We are going to create 2 files that will be used in R in the next step.  

   ##### A. Make blast db on Human FASTA *(only need to run once)*:  
   You will need to set up your dir as follows:   
  * make a fastas subdir - add your human prot fasta to this.  
  * make another subdir called baits - add each of your gene set fastas to this dir.  
 
 Run blastdb script:   
   `./7.A_blastdb.sh`    
   
 output: dir called blastdb  
  
  ##### B. Run BLAST *(for each gene set)*:  
  For each gene set you are investigating, run blast and change the bait file (output of step 6.B.3) for which ever gene set you are on. When blasting - you will only allow 1 hit per query sequence in your bait file.     
   `./7.B_blast.sh ./baits/sensory_percep_light_stim.fa`  
   
   Code in script:   
   `blastp -query $1 -db $fasta -out ${fasta%.}_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 0.00001`   
   
   output: a dir called blastout - within this dir there will be a file called: Homo.fa_ref_blastout
   
  ##### C. Check for any duplicates or missing seqs after blast *(for each gene set)*:  
  Before we check for duplicates in the blastout results, we first need to rename our new blastout directory and blastout file according to the gene set we are working on. Also move the gene_symbol_accid file for this gene set into the blastout dir. For example here is one of my gene set dirs:   
  pwd: hydractinia_transcriptomics/human_seqs_blast/blastout_sens_percep_light_stim  
  Within blastout_sens_percep_light_stim:   
 * Homo.fa_ref_blastout_sens_percep_light_stim   
 * gene_symbol_accid (6.B.3 output file)   
 
   Now that your dir is set up, run 7.C script:    
   `./7.C_check_for_dups_mis.py -a gene_symbol_accid -b Homo.fa_ref_blastout_sens_percep_light_stim`  
   
   This script goes through the blastout file that was just created and compares the accid in the blastout file to the accids in the gene_symbol_accid (6.B.3 output) file to identify if there are any duplicates or missing accids. If duplicates are found, you will need to manually go into the blastout file and remove them - keep the accid with the highest similarity (100) - the script will also tell you how many duplicates are found.   

### 8. Find Hydractinia planula Transcripts in Shared Orthogroups (OGs) with Sensory Genes (In R environment)
  Using the two output files (gene_symbols_accid and blastout) from step 7 and the orthofinder results, we will annotate the Orthogroups with the human gene symbols from the gene sets. We will then identify the Hydractinia planula sequences within the orthogroups of sensory genes for each gene set.   
  
  This step will take place in R so make sure you have the following files on your desktop:   
  * **Rscripts *for each gene set:*** *Hs_Sensory_percep_Light_Stim-updated.Rmd*; *Hs_Sensory_sys_dev-updated.Rmd*; *Hs_TF_activity-updated.Rmd*   
  * **gene_symbol_accid files *for each gene set:*** *gene_symbol_accid_sens_percep_light_stim*; *gene_symbol_accid_sens_sys_dev*; *gene_symbol_accid_total_tf* 
  * **blastout files *for each gene set:*** *Homo.fa_ref_blastout_sens_percep_light_stim*; *Homo.fa_ref_blastout_sens_sys_dev**; *Homo.fa_ref_blastout_total_tf* 
  * **Orthofinder tsv file:** Orthogroups_5-11-20.tsv   
  * **Orthofinder gene count file:** Orthogroups.GeneCount_5-11-20.tsv   

   I have modified portions of an R script made by Sabrina Pankey. **This script has three major parts**, **first** it will read in the orthofinder.tsv file and will re-arrange it so the new dataframe has two columns: seqid and OG (orthogroup) so each row contains 1 seqid with its cooresponding orthogroup number.   
   
   **The second part** of this script identifies the human sequences in the gene set by first reading in the blastout file and the gene_symbol_accid file. It then merges the two files/dataframes by using the function inner_join from dplyr to make a file that contains the gene accid from NCBI, the human seqid from our prot fasta, and the cooresponding orthogroup (OG). From there it then adds another column to this dataframe of the gene symbol.   
   
   **The third part** of this script identifies the Hydractinia planula sequences from the sensory gene orthogroups found in the previous step. A column of the hydractinia seqids are added to this dataframe of the gene accession, human seqid, OG, and gene symbol, one thing to note is that there can be multiple hydractinia seqids in 1 human sensory gene orthogroup or the opposite of no hydractinia seqids in a human sensory gene orthogroup. The script then modifies the headers of the hydractinia seqid, so after TransDecoder our headers looks like this: Hs_planula_t.124..645-1 which was used in orthofinder. In later steps we will need to have the original header so we can use the gene count data from Salmon which is why after Transdecoder we modified our headers. That modification allows us to get back to the original transcript from where the protien model came from. So here the script removes the position in the sequence the protien model came from: Hs_planula_t.124..**645-1** by splitting on the .. so we now have the identity of the original transcript: Hs_planula_t.124. The script then writes out a file to the directory you are working in of this dataframe.    

   What you need to change in this script (in order):   
  * Set your working directory
  * The orthofinder.tsv file name 
  * The blastout file name
  * The gene_symbol_accid file name  
  * The Orthofinder gene count file name
  * The name and path for your output file 
  
   While you are working through this script, I recommend recording key info in a text editor file (one file for each gene set) to make a summary file:  
  * The name of the gene set with the nuber of genes in the gene set
  * The number of human gene models from the second part after making the dataframe of gene acc, seqid, OG, and gene symbol 
  * The number of unique human gene models from the above dataframe
  * The number of Hydractinia gene models found in all of the orthogroups (from part 3)
  * The number of unique hydractinia headers
  * The number of unique Hydractinia OG 

   So in total, you will have 2 key output files:  
  * A summary file made mannually (you will add on to this file in the next step)  
  * The automatically generated file from R that contains the last dataframe:  gene_acc, Homo_seqid, OG, gene_symbol, Hydractinia_seqid   

### 9. Find Significant Differentially Expressed Genes (DEGs) 
  Ultimately, we want to visualize significant gene expression of hydractinia planula sensory genes across developmental days (the 4 days), which we will do in step 10, but first we need to identify significant DEGs by first running EdgeR in R. We will then export the significant pariwise comparisons for each stage into the terminal. Using this, we will compare the headers from the sensory orthorgroups (the output file from step 8) to the significant DEGs from EdgeR to identify signficant actinula genes in the different gene sets. We will then need to run 2 other scripts to re-organize the output info so we can easily add it into the R script of step 10.     

##### A. Run EdgeR script: 9.A_edgeR_hydractinia_4_groups_0.0001.R
   To run EdgeR, we first need to make sure we have our salmon output in our current working directory on our desktop. Make sure all salmon data is in a folder called mapping. You will also need to make a file called libraries_to_stages which will have 2 columns of info, the first column will have the names of all the samples and the second will have the developmental day each sample belongs to. 
  
  This script has 2 major parts. First, it will read in your salmon data (the quant files) and assign them to their corresponding stages based on the libraries_to_stages file. It will then make a multidimensional scaling (MDS) plot based on count data so you can visualize the similarity of each sample in space.  
  
  The second part performs the pairwise comparisons between libraries. For this part you run a function which:   
 * Uses the calcNormFactors function to normalize library sizes (uses trimmed mean of M-values (TMM) method)
 * Estimates Dispersion using a quantile-adjusted conditional maximum likelihood (qCML) method using the estimateCommonDisp funciton.  
 * Examines DEGs at p-value = 0.05 
 
  The Function then generates 4 plots in 1 figure (DGE Exact Test, MDS Plot for Count Data, BCV plot, and Mean-Variance plot) for each comparison. The key output information you need from this analysis to use in the next steps are the files that are generated at each comparison that contain the headers of the significant DEGs with their logFC, logCPM, pvalues, and FDR. 
  
  What you need to change in this script (in order):   
  * Set your working directory
  * Adjust the number of days accordingly for reading in the data, here we have 4 developmental days 
  * Adjust the TagwiseDisp n value: 50/(#samples - #groups) (this is in both functions)
  * The paths and file names for each pairwise comparison   
   
  The output that you will need for the next step: 
 * The output files from each pairwise comparisons - since we have 4 days there will be 6 files generated 
   
##### B. Identify significant DEGs in gene sets (in terminal)  
  Now that we have all of the significant DEGs between each stage, we can now use this to identify which planula genes are significantly expressed in the gene sets (the genes that share orthogroups with human sensory genes). 
  
  In this first part, we are going to run the 9.B_find_sig_degs_in_geneset.py script but first we need to import data and set up our directory in the terminal. I created a new dir called Hydractinia_DEGs within my blastout dir for each of my gene sets from step 7. Import all 6 files of the edgeR output into your dir. I also recommend opening the text summary files that we created in step 8 to add more info. So the 9.B_find_sig_degs_in_geneset.py script takes in 2 files, the output file from step 8 of the headers of interest and one DEG edgeR comparison file. It also takes in the name for the output file. The script will create a dictionary from the step 8 output (since these are the headers that share orthogroups with sensory genes) and then it will compare each header from the edgeR file and will record the edgeR info if there is a matching header. This will then get written to a new file - I recommend copying this info into your summary file. You will have to do this for each edgeR output file (so 6 times) - In the future I want to make it so you only have to run it once (iterate through all the files in the dir using os.walk).   
  
  Run 9.B_find_sig_degs_in_geneset.py *__for each comparison file__*:   
  `./9.B_find_sig_degs_in_geneset.py -a Reduced_Hs_Sensory_Percep_Light_Stim_R_output.txt -b s1vs2_0.0001_4gr -c s1vs2`
   
   Output file: s1vs2_red_DEGs.txt (you will have 6 output files which is why I recommend copying all of the output into your summary file).  
  
  Now we will prep our significant DEG output to be used in step 10. Copy your summary file and import it into the terminal in this dir. We are going to run 2 prep scripts:  
 * The first will go through our summary file and will write out a file with all of our significant DEG headers on one line with quotes around each header seperated by a comma. Input your summary file and the name of the gene set for your output file. Copy the output and add it at the bottom of your summary text document.    
  `./9.C_prep_sig_DEGs_for_heatmap.py -a Summary_DEG_Results_Hs_Sensory_Percep_Light.txt -b sens_percep_light_stim`  
  
 * The second prep script will write out a file that has 2 columns, in the first it has the orthogroup and the second has all of the gene symbols associated with that orthogroup, if there are multiple gene symbols it seperates them with a comma. This script requires 2 input arguments: the output file from step 8 with all of the OG of interest and the gene set name for the output file.   
 `./9.D_get_all_symbols_for_OGs.py -a Reduced_Hs_Sensory_Percep_Light_Stim_R_output.txt -b sens_percep_light_stim`   
  
  
#### 10. Generate Heatmaps of Significant DEGs of Gene Sets (In R environment)
  Now that we have identified the significant actinula DEGs in the gene sets we can visualize their expression by generating 2 heatmaps for each gene set. 
  
  Run the heatmap scripts in R *(run on each gene set)*  
  * 10_Hs_Sensory_percep_Light_Stim_heatmap-updated.R.  
  * 10_Hs_Sensory_Sys_Dev_heatmap-updated.R.  
  * 10_Hs_TF_activitty_heatmap-updated.R  
  
  What you will need for the script:   
 * Salmon output (the mapping dir we used before in edgeR)
 * libraries_to_stages.txt file that we used before in edgeR
 * The R output file from step 8 (contains gene_acc, Homo_seqid, OG, gene symbol, Actinula_seqid)
 * all_symbols_for_OGs output from 9.D_get_all_symbols_for_OGs.py 
 * The output string from 9.C_prep_sig_DEGs_for_heatmap.py (all sig degs in format: "header_1", "header_2", "header_3",...).  
 
 Each script will make 2 heatmaps, the first is the total heatmap with all of the hydractinia genes found in that gene set. The second heatmap contains only significant DEGs from that gene set. This script makes heatmaps using the heatmap.2 function. Each row in the heatmap will be annotated with the hydractinia transcript header and will also have all of the gene symbols assoicated with the orthogroup that header is apart of. There will also only be unique headers displayed in the heatmaps, there are no redundant headers. 
 
 
#### 11. Supplemental Step - Identify Significant Genes Upregulated in Stages 5 and 6   
   In the heatmaps generated in step 10, we found a general pattern that in stages 5 (when the larvae are settled) and 6 (when larvae complete metamorphosis into a juvenile polyp) sensory gene expression is highly downregulated. We want to make sure that this pattern is not an artifact of our analysis, so here we are identifying genes that are significantly upregulated in stages 5 & 6 and identifying their functions by performing a Gene Ontology (GO) analysis using interproscan https://interproscan5-docs.readthedocs.io/en/latest/HowToRun.html and ReviGO http://revigo.irb.hr/.   
   
   ##### A. Find DEGs that are upregulated from EdgeR output   
   First we need to identify all of the headers that are upregulated in both stage 5 and 6. To do this, run the two 11.A scripts which will identify all of the headers that have a positive logFC value from the edgeR output and will write out a file with these headers. Before you run these scripts make a new directory to work in (one for each stage) and copy the edgeR comparison files that involve these 2 stages from step 9: 
  * The stage 5 edgeR files: s1vs5_0.05_6gr; s2vs5_0.05_6gr; s3vs5_0.05_6gr; s4vs5_0.05_6gr   
  * The stage 6 edgeR files: s1vs6_0.05_6gr; s2vs6_0.05_6gr; s3vs6_0.05_6gr; s4vs6_0.05_6gr
  
  Now Run the 2 scripts:  
  `./11.A_find_upregulated_stg_5_genes.py -a s1vs5_0.05_6gr -b s2vs5_0.05_6gr -c s3vs5_0.05_6gr -d s4vs5_0.05_6gr`      
  `./11.A_find_upregulated_stg_6_genes.py -a s1vs6_0.05_6gr -b s2vs6_0.05_6gr -c s3vs6_0.05_6gr -d s4vs6_0.05_6gr` 
  
  Outputs: 
 * stg_5_upregulated_headers.txt  #This file has 263 headers  
 * stg_6_upregulated_headers.txt  #This file has 111 headers   
  
  ##### B. Make FASTA of significant DEGs using selectSeqs
  We need to make FASTA files for both stages so we can run Interproscan in the next step. Run 11.B_selectSeqs.pl script on the stage 5 and 6 outputs from above.   
  `./11.B_selectSeqs.pl -f stg_5_upregulated_headers.txt ../../salmon_3-13-20/actinula_total.ORP.fa-mod.fa >> stg_5.fa` 
  `./11.B_selectSeqs.pl -f stg_6_upregulated_headers.txt ../../salmon_3-13-20/actinula_total.ORP.fa-mod.fa >> stg_6.fa`   
 
 Outputs: 
* stg_5.fa 
* stg_6.fa   

##### C. Run Interproscan, Pull out GO terms, and visualize in ReviGO
 Now that we have our FASTAs of the siginificantly upregulated DEGs for stages 5 & 6, we can run interproscan on both stages to get GO terms. Run interproscan on stage 5 and 6.   
 `sbatch 11.C.1_interpro_stg_5.slurm`    
 `sbatch 11.C.1_interpro_stg_6_slurm`    
 
 Code in scripts:  
 `interproscan -i ./stg_5.fa -d ./interpro_stg_5 -goterms -f TSV`   
 `interproscan -i ./stg_6.fa -d ./interpro_stg_6 -goterms -f TSV` 
 
 Outputs: 
* stg_5.fa.tsv
* stg_6.fa.tsv
 
 Next, we need to pull out the GO terms from the tsv file so we can copy the GO terms into ReviGo.  
 `./11.C.2_make_annotations_file.py -i stg_5.fa.tsv`     
 `./11.C.2_make_annotations_file.py -i stg_6.fa.tsv`   
 
 Outputs: 
* annotations.txt (rename for each stage)  
*For each output file - open it in a text editor. Remove the headers/copy all of the GO terms and paste in a new document - organize them so there is 1 GO term per line - call this new file stg_5_GOs.txt and stg_6_GOs.txt*  
 
 Now that you have all of the GO terms, organized one per line - copy and paste them into ReviGo to visualize functions http://revigo.irb.hr/ 

