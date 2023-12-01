#! /usr/bin/env python3

#This script is one cohesive script (the combo of 1_steps_1-3a_v3.py and 1_steps_3b-4b.py) to prep my raw reads for assembling a reference transcriptome. 
#Our data was sequenced by novogene. We have 6 developmental stages each with 6 replicates
#The goal is to identify the highest quality replicate in each stage to use in the reference transcriptome. We are doing this by using FastQC and calculating basic stats on each
#replicate - which ever rep for each stage gets the highest score will be used in the reference transcriptome (they will be concatenated and used in the R1 and R2 files for assembly).


#import modules        #dont forget to module load linuxbrew/colsa 
import argparse
import subprocess
import os
import gzip
import time
import shutil
from zipfile import ZipFile
import json

#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("--dir", "-d", help="Path to Read directory with sub directories of raw reads - subdir names should be in the format of group_info_additional-info")

args = parser.parse_args()

#time program
t1 = time.time()

#for this data - there are 36 sub dirs, with a total of 72 files (R1 and R2) 

#Step 1: prep files
#1.a) make a fastqc directory and prep for fastqc run (identify groupings of samples) 
#create a directory first then identify groups and all items within that group (the replicates) and write to a file
#for example: {STG1: [STG_1_R1, STG_1_R2 ...] STG2: [STG_2_R1, STG_2_R2...]}

def fastqc_prep (dir):
    group_db = {}
    os.mkdir("fastqc_results")
    for item in os.scandir(dir):
        #print("iterating on item: ", item) #for debug
        if item.is_dir():
            group = item.name.split("_")[0]
            #print("group 1: ", group) #for debug
            group_2 = item.name.split("_")[1]
            #print("group 2: ", group_2) #for debug
            tot_group = group + group_2
            #print(tot_group) #for debug

            group_db.setdefault(tot_group,[])
            group_db[tot_group].append(item.name)
    #print("group db: ", group_db)#for debug

    #Write group structure to an output file
    with open("groups_and_winning_samples.txt", "w") as out_handle:
        out_handle.write("The groups of your samples are as follows: \n ")
        for key in group_db:
            out_handle.write("Group: {0}, samples in that group: {1}\n".format(key, group_db[key]))
    
    return group_db

#call function using the directory given
fastqc_prep_result = fastqc_prep(args.dir)
print("structure of groups and fastqc prep db: ", fastqc_prep_result) #for debug



#Step 2) iterate through subdirectory specified in args.d - copy zipped files and rename them - place in directory specified
try:
    iter_v = 0
    #iterate through the subdirectories in the directory specified - it will copy the files out of the subdirs and change 
    #the names of the read files based on the subdir names 
    for root, dirs, files in os.walk(args.dir): 
        r1_output = "{0}.R1.fastq.gz".format(root)
        #print (r1_output) #debug
        r2_output = "{0}.R2.fastq.gz".format(root)
        #print(r2_output) #debug
  
        iter_v += 1
        print("on iteration: ", iter_v) #for debug
        #In the subdirectory, sort the file, ID if file is R1 or R2 - copy the current file and rename and move it (all while ziped)
        for file in sorted(files): 
            print("on file:",file)

            #id if the current file is R1 or R2 and write to new files 
            if "_1.fq" in file:
                print("On R1 file: ", file)
                #source path and destination path
                cur_d =	os.getcwd()
       	       	#print("cwd: ", cur_d) #debug
                #print("root:", root)#debug
                #print("file: ", file)#debug
                source = "{0}/{1}/{2}".format(cur_d,root, file)
                print("source: ",source)
                destination = "{0}/{1}".format(cur_d, r1_output)
                print("destination: ", destination)
                #copy the file
                dest = shutil.copy(source, destination)
                print("file {0} is copied".format(r1_output))
    

            elif "_2.fq" in file: 
                print("On R2 file: ", file)
                #source path and destination path
                cur_d =	os.getcwd()
       	       	#print("cwd: ", cur_d) #debug
                #print("root: ", root) #debug
                #print("file", file) #debug
                source = "{0}/{1}/{2}".format(cur_d,root, file)
                print("source: ",source)
                destination = "{0}/{1}".format(cur_d, r2_output)
                print("destination; ", destination)
                #copy the file
                dest = shutil.copy(source, destination)
                print("file{0} is copied".format(r2_output))
                

except IOError as err:
        print("problem reading or writing/appending file:", err)
        
print("heading to step 3 - fastqc")

#Step 3 run fastqc and calculate scores  
#3a) run fastqc using a function - this will run in the directory specified (args.d) using the zipped 
#files copied from step 2 - output will be in the fastqc dir which is outside the dir specified 
def run_fastqc():
    os.chdir(args.dir)
    #print(os.getcwd()) #for debug
    for item in os.scandir("."):
        if item.is_file(): 
            if ".R1.fastq.gz" in item.name: 
                #print("item: ", item)#debug
                #print("item.name: ", item.name)#debug
                print("Running fastqc on: {0}".format(item.name))
                result = subprocess.run("fastqc {0} -o ../fastqc_results".format(item.name), shell=True)
                print("fastqc complete on {0}".format(item.name))
                
            elif ".R2.fastq.gz" in item.name: 
                print("Running fastqc on: {0}".format(item.name))
                result = subprocess.run("fastqc {0} -o ../fastqc_results".format(item.name), shell=True)
                print("fastqc complete on {0}".format(item.name))


    print("All fastqc is complete")

#call run_fastqc function
result_1 = run_fastqc()

print("Steps 1-3a (fastqc) is complete - move on to next steps")



#3b) iterate through keys and values of group_db to get the base name of files to use in a path to open fastqc file and collect info    
try: 
    def stats(db):
        #make empty dictionaries to keep track of stats values with keys 
        max_db = []
        min_db = []
        stats_db = {}  #order of columns in db(they are in a list): total seqs, Seqs length, quality (avg of mean) 
        #each key will be 1 file ex: {STG_1_R1.R1 : [total seqs, seqs length, quality(avg of mean) ], STG_1_R1.R2 : [total seqs, seqs length, quality ] ...}
    
        os.chdir("../fastqc_results")
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
results = stats(fastqc_prep_result) #if running in 1 full script use this call 
#results = stats(group_db) #if doing in steps use the group_db variable 


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
                #print("key in calcs db: ", sc_key)
                if sc_key.startswith(value) and sc_key.endswith(".R1"): 
                    score_list_r1.append(score_db[sc_key])
                    name_list_r1.append(sc_key)
                  
                    print("on value: ", value)
                    print("score list: ", score_list_r1)
                    print("name list: ", name_list_r1)

        #id the highest score per group using the counter variable, the index variable keeps track of the position
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
winners = select_best_reads(fastqc_prep_result, results)
print("returned winning list from function: ", winners)


#make a list of the coresponding R2  winning reads 
winners_2 = []
for name in winners: 
    new_name = name.replace(".R1", ".R2")
    winners_2.append(new_name)
print("winners 2 ", winners_2)

#returned collection from step 4a: 
#winners = ['STG_3_R2.R1', 'STG_5_R2.R1', 'STG_2_R4.R1', 'STG_6_R4.R1', 'STG_1_R6.R1', 'STG_4_R3.R1']
#winners_2 = ['STG_3_R2.R2', 'STG_5_R2.R2', 'STG_2_R4.R2', 'STG_6_R4.R2', 'STG_1_R6.R2', 'STG_4_R3.R2']


#Write winners to be concatenated to output file
with open("../groups_and_winning_samples.txt", "a") as out_handle: 
    out_handle.write("\nThe winning (highest quality = representitive) samples from each group that will be concatenated in your total_R1 and total_R2 files are: \n")
    for sample in winners: 
        out_handle.write(sample + "\t")
    
    out_handle.write("\n")    
        
    for sample in winners_2: 
        out_handle.write(sample + "\t")   


#Step 4 b) Concatenate the winning reads into Total R1 and Total R2 files (zipped) 
def make_totals_files(samp_dir, win_list, win_list2):
    #change into the correct dir, open new files and concatenate the highest scoring reps for each stage in order 
    #os.chdir("./fastqc_results") #only for testing - this is where we are after the previous code 
    print("Beginning concatenation of Total R1 and Total R2 files, cwd: ", os.getcwd()) #for debug
    os.chdir("../{0}".format(samp_dir))
    print("cwd after cd: ", os.getcwd()) #for debug
    
    #open two zipped files to write to - total_R1 and total_R2
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
                    
    print("total R1 and R2 files have been made!!")

#Call function 
make_totals_files(args.dir, winners, winners_2)

#time program
t2 = time.time()
secs_time = t2-t1
int(secs_time)
tot_time = secs_time/60

print("Total running time: {0} minutes ".format(tot_time)) #in minutes
