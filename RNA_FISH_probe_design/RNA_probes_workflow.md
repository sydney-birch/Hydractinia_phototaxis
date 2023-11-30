# RNA FISH Probe Design Workflow

This file is the step-by-step instructions for designing RNA Fluorescent in situ hybridization (FISH) probes for actinulae larvae (Hydrozoa) using Stellaris Custom RNA FISH Probe sets https://www.biosearchtech.com/products/rna-fish/custom-stellaris-probe-sets. I designed these probes based on the results of the gene expression work and outputs which will be refered to here. Scripts for this section can be found in the directory: Actinula_Paper/RNA_FISH_probe_design/Probe_design_scripts

## Process overview:
1. Open the R output files for your gene set of interest (generated at step 8 of the Transcriptome analysis)
2. Make a list of the OGs of insterest and their gene symbols and actinula sequence IDs
3. Open the Rmd script - go through until you get the full hydractinia seq header with seq location - add to file #go to line 146
4. Use script to open salmon data and get expression values --> makes a csv
5. Open csv and select the highest expressed headers for each OG
6. Open alignment files for OGs of interest - blastp all seqs and compare to expression data - select best seqs
7. Pull out all FASTA seqs, get exact seq (use snapgene) & blastx
8. Choose best & make probes in stallaris


### 1. Open R.md output file from step 8 of the transcriptome workflow
For the geneset of interest, open the output file from the R.md script generated in step 8 of the transcriptome workflow:

- Sensory Perception of light stimulus: 8_Reduced_Hs_Sensory_Percep_Light_Stim_R_output_genes.txt

Each file contains the following information for all hydractinia sequences that are found in any of our annotated OGs for the cooresponding geneset: gene_acc, Homo_seqid, OG, gene symbol, Actinula_seqid

### 2. Make a list of the OGs of interest
Next, make a list of all OGs you are interested in making probes for

RRH,OPN5,RHO,OPN4 OG0000063
CRX;VSX1;VSX2;RAX2  OG0000012
CNGA3;CNGB3	OG0001196

### 3. Get full sequence headers with sequence locations
A. Open the Rmd scripts for the geneset of interest and go through the script until line 146.
B. Open the actinula_acc_apo_OGs_acc_sym dataframe and search the OG of interest
C. For each OG, write out the full hydractinia sequence header which contains the sequence location - we will need that later in Step 7.

Hs_planula_t.11713..435-1	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.14852..1060-236	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.20655..3-353	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.21263..1204-488	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.21264..255-572	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.22812..451-65	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.23688..188-1144	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.24723..1175-555	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.24723..9401-9096	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.3778..394-2	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.40086..1424-483	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.40948..792-103	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.47300..1369-698	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.50049..190-774	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.74880..126-824	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.74881..1377-721	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.85107..1378-566	CRX;VSX1;VSX2;RAX2	OG0000012 
Hs_planula_t.87239..1008-103	CRX;VSX1;VSX2;RAX2	OG0000012 

Hs_planula_t.10798..710-3	RRH;OPN5;RHO;OPN4	OG0000063
Hs_planula_t.29380..1208-315	RRH;OPN5;RHO;OPN4	OG0000063
Hs_planula_t.75305..60-938	RRH;OPN5;RHO;OPN4	OG0000063
Hs_planula_t.88569..1379-489	RRH;OPN5;RHO;OPN4	OG0000063

Hs_planula_t.15325..1-420	CNGA3;CNGB3	OG0001196
Hs_planula_t.35573..756-3677	CNGA3;CNGB3	OG0001196


### 4. Make a CSV of gene expression values for sequences of interest
In the terminal, make a probes directory where you can easily navigate to your salmon dir from the transcriptome analysis.
`mkdir hydractinia_transcriptomics/hydractinia_probes`

Next, in the hydractinia_probes dir, make dirs for each of your OGs of interest:

`mkdir hydractinia_transcriptomics/hydractinia_probes/opsins ` 
`mkdir hydractinia_transcriptomics/hydractinia_probes/CNGs`  
`mkdir hydractinia_transcriptomics/hydractinia_probes/CRXs `

In each dir, make a headers of interest file and copy the information we listed out in Step 3. Each line should contain the following information in order for however many headers you have in each OG: complete header (described in the above step), the gene symbol, and OG. (Each item should be tab seperated)

`nano CRXs_headers_of_interest_OG0000012 `
`nano OPNs_headers_of_interest_OG0000063`
`nano CNGA3_headers_of_interest_OG0001196`

Next, Run the make_csv_from_OF_and_quant_files_and_sort_TPM.py on each OG of interest. This script will navigate to the salmon dir and will generate a CSV file that has all of the count information for each of the headers of interest, and will sort the file to have the highest expressed transcripts at the top. This script requires 4 pieces of info at execution: 1) the headers of interest file, 2) the path to the salmon dir with your read quantification info, 3) the first part of each group of transcripts (STG_1_R1, STG_2_R2 - the program can ignore the STG designation and use the important piece of info), and 4) the name of the output files.

Opsins:
Run script
`./make_csv_from_OF_and_quant_files_and_sort_TPM.py -a OPNs_headers_of_interest_OG0000063 -b ../../salmon/ -c STG_ -d OPNs_OF_expression_info`

output:
OPNs_OF_expression_info.csv
sorted_TPM_file.csv #rename: OPNs_hydractinia_sorted_TPM_file.csv

            #copy to desktop  
	`scp sjb1061@premise.sr.unh.edu:./hydractinia_transcriptomics/Probe_design_5-6-21/OPNs/OPNs_hydractinia_sorted_TPM_file.csv ./`  
 
CRXs:
Run script
`./make_csv_from_OF_and_quant_files_and_sort_TPM.py -a CRXs_headers_of_interest_OG0000012  -b ../../salmon/ -c STG_ -d CRXs_OF_expression_info`

output:
CRXs_OF_expression_info.csv
sorted_TPM_file.csv #rename CRXs_hydractinia_sorted_TPM_file.csv

    #copy to desktop  
	`scp sjb1061@premise.sr.unh.edu:./hydractinia_transcriptomics/Probe_design_5-6-21/CRXs/CRXs_hydractinia_sorted_TPM_file.csv ./`  


Run on CNGs:
Run script
`./make_csv_from_OF_and_quant_files_and_sort_TPM.py -a CNGA3_headers_of_interest_OG0001196 -b ../../salmon/ -c STG_ -d CNGA3_OF_expression_info`

output:
CNGA3_OF_expression_info.csv
sorted_TPM_file.csv #rename CNGA3_hydractinia_sorted_TPM_file.csv

	#copy to desktop  
	`scp sjb1061@premise.sr.unh.edu:./hydractinia_transcriptomics/Probe_design_5-6-21/CNGA3/CNGA3_hydractinia_sorted_TPM_file.csv ./`   


### 5. Select the highest expressed actinula sequences for each OG
Next we will select the highest expressed hydractinia transcripts to investigate as our probe seequences.

1. For each of your OGs of interest, open the hydractinia_sorted_TPM_file.csv
2. Next, filter the table by the first unique transcript header and record the Sample stage/group with the highest TPMs (Transcript per million). For hydractinia, stages 3 and 4 are when larvae are phototactice. We want to make probes out of the highest expressed transcripts during stages 3 and 4. Do this filter step for the top 2-4 highhest expressed transcripts and recrod the unique transcript header, and the TPM range of stages 3 and 4.

### 6. Blastp highly expressed sequences
Next, we will look at the alignments of the sequences in each OG, which can be found in the OrthoFinder MultipleSequenceAlignments output dir. Create a new dir that will contain all of your alignments outside of the MultipleSequenceAlignments dir (side note: don't ls in MultipleSequenceAlignments). Then transfer alignment to desktop to view alignment

   `mkdir OrthoFinder/Results/OGs_of_interest_alignments`

   ```
   cp MultipleSequenceAlignments/OG0000063.fa OGs_of_interest_alignments_3-18-21 #OPN
   cp MultipleSequenceAlignments/OG0001196.fa OGs_of_interest_alignments_3-18-21 #CNG
   cp MultipleSequenceAlignments/OG0000012.fa OGs_of_interest_alignments_3-18-21 #CRX
   ```
 
Next, open the file in seaviewer (or other alignment viewer) and compare the sequences - make note of potential good/bad seqs based on alignment. Then copy the hydractinia sequences into a text editor and remove all dashes (find/replace dash with nothing). Now blastp the sequences in NCBI and record hits and potential sequences for probes. This is a reciprocal blast to try to verify highly expressed seqences to make probes.
```
#CRXs:
*Hs_planula_t.47300..1369-698	STG_1	18.77 - 11.6	Select seq gb|ALJ33542.1|	Rx [Clytia hemisphaerica]
Hs_planula_t.85107..1378-566	STG_1	10.9 - 2.6		Select seq ref|NP_001302439.1|	retinal homeobox protein Rx-B-like [Hydra vulgaris]
Hs_planula_t.50049..190-774	STG_2	10.9-2.9		Select seq ref|XP_002168027.1|	PREDICTED: aristaless-related homeobox protein-like [Hydra vulgaris]	

#OPNs:
*Hs_planula_t.88569..1379-489	STG_2	27.68 - 13.7	Select seq gb|AGB67492.1|	c-like opsin [Tripedalia cystophora]
Hs_planula_t.29380..1208-315	STG_2	21.3 - 8.5		Select seq dbj|BAG80696.1|	opsin [Carybdea rastonii] and Select seq gb|AGB67492.1|	c-like opsin [Tripedalia cystophora]
Hs_planula_t.75305..60-938	STG_4	10.6 - 5.1		Select seq dbj|BAF95843.1|	opsin [Podocoryna carnea] and Select seq ref|XP_012555438.1|	PREDICTED: melanopsin-B-like [Hydra vulgaris]

#CNGs:
*Hs_planula_t.35573..756-3677	STG_3	2.4 - 0.69		Select seq ref|XP_012555740.1|	PREDICTED: cyclic nucleotide-gated cation channel alpha-3-like isoform X1 [Hydra vulgaris]
Hs_planula_t.15325..1-420	STG_2	0.9 - 0.0		Select seq ref|XP_012555744.1|	PREDICTED: cyclic nucleotide-gated cation channel alpha-3-like isoform X2 [Hydra vulgaris]
```

### 7. Pull out FASTA sequences for all transcripts of interest and get exact sequence locaiton of transcript
Next, we will pull out the sequences across all OGs that we are interested in making probes for.

First, make a header file of all the headers you want to pull sequences out call it: nano actinula_headers_for_probes
Copy all headers of interest into the hydractinia_expression_headers  file:

#hydractinia_expression_headers  file should look like this: 
#CRXs
Hs_planula_t.47300..1369-698
Hs_planula_t.85107..1378-566
#OPNs	
Hs_planula_t.88569..1379-489
Hs_planula_t.29380..1208-315
#CNGA3
Hs_planula_t.35573..756-3677

Next, run selectSeqs.pl to get all of the fasta nuc seqs for these headers. You will give this script the headers file, the transcriptome fasta, and the output file name:

	./selectSeqs.pl -f hydractinia_expression_headers ./hydractinia_total.ORP.fasta-mod.fa >> hydractinia_highest_expression_sensory_seqs.fa
 
Save this fasta file on your local computer and open it with a sequence viewer, here I use the free version of SnapGene. In SnapGene, each header will have its own page, go through each one, click on the translation option which shows the reading frame, and copy the specified sequence region from the location portion of the protien header. This will be the exact sequence you use to make probes.

- Take note of the bp and GC content in your notes file. You will do this for each of your sequences.

### 8. Design probes in stellaris and order
Now that you have a fasta file of nucleotide sequences of your exact region for probes, we are going to make a custom probeset using Stellaris: https://www.biosearchtech.com/products/rna-fish/custom-stellaris-probe-sets. I recommend reading over the design information provided by stellaris.

Select the Stellaris rna fish probe designer option and start design.
Enter in the name of your probe set name, for organism select other, masking level change to 2, Max number of probes keep at 48. For the first test keep oligo length at 20 and min spacing length (nt) at 2 - alter these according to the recomendations depending on the number of probes generated. Next, paste in your target sequence and click design probes.
This will generate a probe set. Stellaris recommends a minimum of 25 oligos for a single probe set. If the count number is too low, alter the length and min spacing length according to their recomendations for troubleshooting/designing: https://blog.biosearchtech.com/considerations-for-optimizing-stellaris-rna-fish-probe-design.
Next, select your stellaris dye and order! Check out https://www.biosearchtech.com/support/education/stellaris-rna-fish/dyes-and-modifications-for-stellaris.
You wll go through this with all of your potential sequences and choose the best probe set to order for each OG, we chose:

>Ec_actinula_t.72976..146-1048(+)_opsin_C		903bp 	48% GC	actinula_opsin_C; 34 probes; C3-Fluorescein; 5nmol 
>Ec_actinula_t.17544..7629-1057_(-)_Piezo		6573bp 	46% GC	Actinula_Piezo; 48 probes; Quasar 670; 5 nmol total: 782.20 (NOTE: I Used the - strand (says to use sense strand))
>Ec_actinula_t.81080..203-2761_(+)_PKD2L1_A		2559bp 	48% GC	Actinula_PKD2L1_A; 48 probes; Tamara c9; 5nmol 
>Ec_actinula_t.66208..2940-1282_(-)_PKD1L3_A 		1659bp 	45% GC	Actinula_PKD1L3_A; 48 probes; Tamara C9; 5nol (NOTE: I used the - strand )
>Ec_actinula_t.66269..63-3458_(+)_TRPA_A		3396bp 	46% GC	Actinula_TRPA_A; 48 probes; Quasar 670; 5nmol 

You have now completed the RNA FISH Probe Design!!
