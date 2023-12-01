#! /usr/bin/env python3
#This script is the first of 2 scripts to prep my raw reads for assembling a reference transcriptome. Our data was sequenced by novogene. We have 6 developmental stages each with 6 replicates
#The goal is to identify the highest quality replicate in each stage to use in the reference transcriptome. We are doing this by using FastQC and calculating basic stats on each
#replicate - which ever rep for each stage gets the highest score will be used in the reference transcriptome (they will be concatenated and used in the R1 and R2 files for assembly).

#import modules 
import argparse
import subprocess
import os
import gzip
import time
import shutil

#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("--dir", "-d", help="Path to directory with sub directories of raw reads - subdir names should be in the format of group_additional-info")

args = parser.parse_args()


#Step 1: prep files
#1.a) make fastqc directory and prep for fastqc run - create a dictionary with a group and all items within that group (the replicates)
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
    return group_db

#call funciton
fastqc_prep_result = fastqc_prep(args.dir)
print("structure of fastqc prep db: ", fastqc_prep_result) #for debug



#Step 2) iterate through subdirectories to concatenate raw reads(without unziping) - (don't need to fix R2 headers here)
try:
    iter_v = 0
    #open 2 output files for writing so you are iterating through directories only once - helps performance 
    for root, dirs, files in os.walk(args.dir): 
        r1_output = "{0}.R1.fastq.gz".format(root)
        #print (r1_output) #debug
        r2_output = "{0}.R2.fastq.gz".format(root)
        #print(r2_output) #debug
  
        iter_v += 1
        print("on iteration: ", iter_v) #for debug
        for file in sorted(files): 
            print("on file:",file)

                #id if the current file is R1 or R2 and write to new files 
            if "_1.fq" in file:
                print("On R1 file: ", file)
                #source path and destination path
                cur_d =	os.getcwd()
       	       	#print("cwd: ", cur_d)
                #print("root:", root)
                #print("file: ", file)
                source = "{0}/{1}/{2}".format(cur_d,root, file)
                print("source: ",source)
                destination = "{0}/{1}".format(cur_d, r1_output)
                print("destination: ", destination)
                #copy
                dest = shutil.copy(source, destination)
                print("file {0} is copied".format(r1_output))
    

            elif "_2.fq" in file: 
                print("On R2 file: ", file)
                #source path and destination path
                cur_d =	os.getcwd()
       	       	#print("cwd: ", cur_d)
                #print("root: ", root)
                #print("file", file)
                source = "{0}/{1}/{2}".format(cur_d,root, file)
                print("source: ",source)
                destination = "{0}/{1}".format(cur_d, r2_output)
                print("destination; ", destination)
                #copy
                dest = shutil.copy(source, destination)
                print("file{0} is copied".format(r2_output))
                

except IOError as err:
        print("problem reading or writing/appending file:", err)
        
print("heading to step 3 - fastqc")

#Step 3 run fastqc and calculate scores  
#3a) run fastqc using a function   
def run_fastqc():
    os.chdir(args.dir)
    #print(os.getcwd()) #for debug
    for item in os.scandir("."):
        if item.is_file(): 
            if ".R1.fastq.gz" in item.name: 
                #print("item: ", item)
                #print("item.name: ", item.name)
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
