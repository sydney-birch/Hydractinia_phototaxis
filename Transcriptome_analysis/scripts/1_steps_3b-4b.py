#! /usr/bin/env python3

#The previous script organized group/rep structure and ran FastQC, now this script will calculate stats using the FastQC output to determine which rep for each stage
#is the highest quality rep to be used in the reference transcriptome. For this script I did have to mannually hard code the collection returend from step 1 and assigned it 
#to the variable group_db on line 23. I will add a script that is the full process of all 3 scripts that doesn't require the user to hard code a collection. 

#import modules        #dont forget to module load linuxbrew/colsa 
import argparse
import subprocess
import os
import gzip
import time
import shutil
from zipfile import ZipFile

#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("--dir", "-d", help="Path to directory with sub directories of raw reads - subdir names should be in the format of group_additional-info")

args = parser.parse_args()

#returned collection from step 1:  
group_db = {'STG3': ['STG_3_R6', 'STG_3_R3', 'STG_3_R4', 'STG_3_R5', 'STG_3_R1', 'STG_3_R2'], 'STG5': ['STG_5_R5', 'STG_5_R1', 'STG_5_R2', 'STG_5_R3', 'STG_5_R6', 'STG_5_R4'], 'STG2': ['STG_2_R4', 'STG_2_R2', 'STG_2_R1', 'STG_2_R6', 'STG_2_R5', 'STG_2_R3'], 'STG6': ['STG_6_R6', 'STG_6_R1', 'STG_6_R2', 'STG_6_R5', 'STG_6_R4', 'STG_6_R3'], 'STG1': ['STG_1_R5', 'STG_1_R2', 'STG_1_R1', 'STG_1_R6', 'STG_1_R3', 'STG_1_R4'], 'STG4': ['STG_4_R6', 'STG_4_R5', 'STG_4_R2', 'STG_4_R3', 'STG_4_R4', 'STG_4_R1']}



#3b) iterate through keys and values of group_db to get the base name of files to use in a path to open FastQ files    
try: 
    def stats(db):
        #make empty dictionaries to keep track of stats values with keys 
        max_db = []
        min_db = []
        stats_db = {}  #order of columns in db(they are in a list): total seqs, Seqs length, quality (avg of mean) 
        #each key will be 1 file ex: {STG_1_R1.R1 : [total seqs, seqs length, quality(avg of mean) ], STG_1_R1.R2 : [total seqs, seqs length, quality ] ...}
    
        os.chdir("fastqc_results")
        print("cwd: ", os.getcwd())
    
        #iterate through the base names - open up the R1 and R2 using this base name and add 
        #total seq, seq length and calc avg of mean for the quality and add to lists and dbs also record min and max for each file in lists
        for key in db:
            for value in db[key]:
                print ("on value: ", value) #for debug
                #print(os.getcwd()) #for debug
   
                with ZipFile("./{0}.R1_fastqc.zip".format(value), "r") as in_handel_1:
                    with in_handel_1.open("{0}.R1_fastqc/fastqc_data.txt".format(value), "r") as internal_handel:
                        #cur_file = in_handel_1.read("{0}.R1_fastqc/fastqc_data.txt".format(value))
                        print("On file {0}.R1_fastqc".format(value))
                        avg_list = []
                        trigger = False
          		
                        for line in internal_handel:
                            line = line.decode("ascii")
                            #print("line: ", line) 
                            line = line.rstrip()
                            sp_line = line.split("\t")
                            #print("split line: ", sp_line) #for debug
		
                            if sp_line[0].startswith("Total Sequences"):
                                stats_db.setdefault("{0}.R1".format(value), [])
                                stats_db["{0}.R1".format(value)].append(sp_line[1])
                                #print("Tot seq sp line: ", sp_line[1]) # for debug
		
                            if sp_line[0].startswith("Sequence length"):
                                stats_db.setdefault("{0}.R1".format(value), [])
                                stats_db["{0}.R1".format(value)].append(sp_line[1])
                                #print("sq len sp line ", sp_line[1]) #for debug
		
                            #Extract the mean and take the average of the mean and add to the stats_db
                            if sp_line[0] == "#Base":
                                trigger = True
		
                            elif trigger:
                                if sp_line[0] == ">>END_MODULE":
                        	    #calculate the average of the means and append to list then find the lowest min and highest max for all files 
                                    avg_of_mean = sum(avg_list)/len(avg_list)
                                    #print("avg:", avg_of_mean) #for debug
                                    stats_db["{0}.R1".format(value)].append(avg_of_mean)
		
                                    current_min = min(avg_list)
                                    min_db.append(current_min)
                                    #print("min", current_min) #for debug
		
                                    current_max = max(avg_list)
                                    max_db.append(current_max)
                                    #print("max", current_max) #for debug
		
                                    break
		
                                else:
                                    avg_list.append(float(sp_line[1]))
		
                            #print("avg list for iter: ", avg_list) #for debug
                        #print("max db at end of iteration: ", max_db) #for debug
                        #print("min db at end of iteration: ", min_db) #for debug
		
                #Iterate through R2 file
                with ZipFile("./{0}.R2_fastqc.zip".format(value), "r") as in_handel_2:
                    with in_handel_2.open("{0}.R2_fastqc/fastqc_data.txt".format(value), "r") as internal_handel_2:
                        print("On file {0}.R2_fastqc".format(value))
                        #cur_file = in_handel_1.read("{0}.R1_fastqc/fastqc_data.txt".format(value))
                        avg_list = []
                        trigger = False
          		
                        for line in internal_handel_2:
                            line = line.decode("ascii")
                            #print("line: ", line) 
                            line = line.rstrip()
                            sp_line = line.split("\t")
                            #print("split line: ", sp_line) #for debug
		
                            if sp_line[0].startswith("Total Sequences"):
                                stats_db.setdefault("{0}.R2".format(value), [])
                                stats_db["{0}.R2".format(value)].append(sp_line[1])
                                #print("Tot seq sp line: ", sp_line[1]) # for debug
		
                            if sp_line[0].startswith("Sequence length"):
                                stats_db.setdefault("{0}.R2".format(value), [])
                                stats_db["{0}.R2".format(value)].append(sp_line[1])
                                #print("sq len sp line ", sp_line[1]) #for debug
		
                	    #Extract the mean and take the average of the mean and add to the stats_dbb
                            if sp_line[0] == "#Base":
                                trigger = True
		
                            elif trigger:
                                if sp_line[0] == ">>END_MODULE":
                        	    #calculate the average of the means and append to list then find the lowest min and highest max for all files 
                                    avg_of_mean = sum(avg_list)/len(avg_list)
                                    #print("avg:", avg_of_mean) #for debug
                                    stats_db["{0}.R2".format(value)].append(avg_of_mean)
		
                                    current_min = min(avg_list)
                                    min_db.append(current_min)
                                    #print("min", current_min) #for debug
		
                                    current_max = max(avg_list)
                                    max_db.append(current_max)
                                    #print("max", current_max) #for debug
		
                                    break
		
                                else:
                                    avg_list.append(float(sp_line[1]))
		
                            #print("avg list for iter: ", avg_list) #for debug
                        #print("max db at end of iteration: ", max_db) #for debug
                        #print("min db at end of iteration: ", min_db) #for debug
		

        print("stats db at end of value collection: ", stats_db)
        #print("max db at end: ", max_db) #for debug
        #print("min db at end: ", min_db) #for debug


#Step3 c) perform calculations to get the final score 
        total_min = min(min_db)
        total_max = max(max_db)
        #print("Total min across all files: ", total_min) #for debug
        #print("Total max across all files: ", total_max) #for debug
    
        #use the stats_db to calculate the final score - add that final score to a db with the file name as the key
        calcs_db = {}
        #each key will be 1 file ex: {STG_1_R1.R1 : [total seqs, seqs length, quality(avg of mean) ], STG_1_R1.R2 : [total seqs, seqs length, quality ] ...}
        for key in stats_db:
            #print("key in stats db: ", key)
            calcs_db.setdefault(key, 0)
            num_seqs = stats_db[key][0]
            int(num_seqs)
            #print("Number of seqs: ", num_seqs)
            read_length = stats_db[key][1]
            int(read_length)
            #print("Read length: ", read_length)
            quality = stats_db[key][2]
            #print("Read quality: ", quality)

            num_nucs = int(num_seqs) * int(read_length)
            #print("Number of nucs: ", num_nucs) #for debug
            normalized_qual = (quality-total_min)/(total_max-total_min)
            #print("Normalized quality: ", normalized_qual) #for debug
            final_score = num_nucs*normalized_qual
            #print("Final Score: ", final_score) #for debug
            calcs_db[key]+=final_score

        print("final score - calcs db: ", calcs_db) #for debug

        return calcs_db




except IOError as err:
    print("problem reading or writing/appending file:", err)


#Call stats function 
#results = stats(fastqc_prep_result) #if running in 1 full script use this call 
results = stats(group_db) #if doing in steps use the group_db variable - need to hardcode the db in

print("Going to Step 4 a")

#Step 4 a) identify the best reads to use for reference transcriptome 
def select_best_reads(group_db, score_db):
    winners_r1 = []
    #keep track of the name and the score by using 2 lists - names and scores will be in the same positions - ordered 
    for key in group_db:
        score_list_r1 = [] #score 
        name_list_r1 = [] #cooresponding name of file

        for value in group_db[key]:  
            #print("group db key name: ", key)
            #print("group db value: ", value)

            for sc_key in score_db: 
                print("key in calcs db: ", sc_key)
                if sc_key.startswith(value) and sc_key.endswith(".R1"): 
                    score_list_r1.append(score_db[sc_key])
                    name_list_r1.append(sc_key)
                  
                    print("on value: ", value)
                    print("score list: ", score_list_r1)
                    print("name list: ", name_list_r1)

        #id the highest score per group using the counter variable, the index variable keeps track of where the position
        #of the highest score - once finished iterating through lists it will append the highest scoring name to the winners list per group 
        counter = 0
        idx = 0
        for item in range(len(score_list_r1)):
            if counter < score_list_r1[item]:
                counter = score_list_r1[item]
                print("counter", counter)
                idx = item
                print("idx", idx)
            else:
                continue


        winner_r1_name = name_list_r1[idx] 
        winners_r1.append(winner_r1_name)
        print("winners list", winners_r1)  
        
        print("winner names:", winner_r1_name)
    
    return winners_r1

#Call function
winners = select_best_reads(group_db, results)
print("returned list from function: ", winners)


#make a list of the coresponding R2  winning reads 
winners_2 = []
for name in winners: 
    new_name = name.replace(".R1", ".R2")
    winners_2.append(new_name)
print("winners 2 ", winners_2)


print("Going to step 4b")
#Step 4 b) Concatenate the winning reads into Total R1 and Total R2 files (zipped)
def make_totals_files(samp_dir, win_list, win_list2):
    #change into the correct dir, open new files and concatenate the highest scoring reps for each stage in order  
    #os.chdir("./fastqc_results") #only for testing - this is where we are after the previous code 
    print("Beginning concatenation of Total R1 and Total R2 files, cwd: ", os.getcwd()) #for debug
    os.chdir("../{0}".format(samp_dir))
    print(os.getcwd()) #for debug
    
    #Open two zipped files to write to - total_R1 and total_R2 
    with gzip.open("total_R1.fastq.gz", "wb") as out_handel_1:
        with gzip.open("total_R2.fastq.gz", "wb") as out_handel_2:

            for idx in range(len(win_list)): 
                match_name_1 = win_list[idx]
                print("idx: ", idx)
                print("match name 1: ", match_name_1)
                match_name_2 = win_list2[idx]
                print("match name 2: ", match_name_2)

                #write total r1 file - read in each line and write out each line
                with gzip.open("{0}.fastq.gz".format(match_name_1), "rb") as in_handel_1:
                    print("writing for file:", match_name_1)
                    for line in in_handel_1:
                        out_handel_1.write(line)
                        
                #write total r2 file - read in each line and write out each line
                with gzip.open("{0}.fastq.gz".format(match_name_2), "rb") as in_handel_2:
                    print("writing for file", match_name_2)
                    for line in in_handel_2:
                        out_handel_2.write(line)
                    
    print("total R1 and R2 files have been made")

#Call function 
make_totals_files(args.dir, winners, winners_2)
