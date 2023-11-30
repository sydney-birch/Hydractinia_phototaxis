#! /usr/bin/env python3

#Don't forget to: 
#module load anaconda/colsa
#conda find python3 -a 

#This script is the modified version of 5.C_make_csv_of_blastout_and_quant_files_and_sort_TPM.py 
#to use after orthofinder when you are trying to find highly expressed transcripts 
#from a file of headers of interest using your salmon quant.sf file. The script is set up as functions
#and goes through two parts for each sequence - 1) it will write out the input of interest from the input file
#2)Then using that info it will use another function to go to the salmon data and basically run the extract quant
#file py script and write that out to the same file 
#it uses these 2 main functions for each line in the headers of interest file - make another script after this one
#that will identify the highest expressed transcript for each header of interest from each sample used in salmon

#for the input file - make a tab delimited file with the full header (w/seq location), then gene symbol, then OG
#make the file after going through the Rmd script where it makes a dataframe with all the actinula seqs, the OG, and gene symbols

#import modules
import argparse
import os
import pandas as pd
import csv
import operator

#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="input file: tab deliminated file of headers, gene symbol, and OG of interest")
parser.add_argument("-b", help="input path to dir: write the path to the dir that contains all of your salmon output - no quotations (ex: ../../../3_salmon/)")
parser.add_argument("-c", help="groupings of samples in salmon dir (all samples start with HY_ or another example all samples start with STG_)")
parser.add_argument("-d", help="name of output file")

args = parser.parse_args()


#Prep: open file for writing and write the format of the output file
with open("{0}.csv".format(args.d), "w", newline="") as out_handle:
    heading = "Unique_Transcript_Header,Orthogroup,Gene_Symbols,Transcript_Header,Seq_Location_of_Transcript,Sample,TPM,NumReads"
    #split string into list to write to csv
    split_heading = heading.split(",")
    print("split headings list: ", split_heading)
    writer = csv.writer(out_handle, delimiter=",")
    writer.writerow(split_heading)



#Part 1: define function for salmon info
def get_salmon_info(header,dir_path):
    #change dir to your salmon data 
    cwd = os.getcwd()
    #print("current working dir before change: ", cwd)
    
    os.chdir(dir_path)
    #print("current working dir after change: ", os.getcwd())
    
    #This block goes through the directory given and sets the dirs for each sample as keys 
    #in a dictionary and makes a list for the value 
    try:
        counts_db = {}
        for item in os.scandir("."): 
            if item.is_dir(): 
                if args.c in item.name: 
                    #print("current dir name: ", item.name)

                    counts_db.setdefault(item.name, [])
        #print("final db: ", counts_db)

    except IOError as err:
        print("problem reading or writing/appending file:", err)     

    #This section will iterate through the keys of the dictionary (the sample dirs from salmon) - for each key it will open 
    #The quant.sf files and search for the current header, it will record the data in a list 
    #and will write out the info to the output file (with the info from the blastout file) 
    try: 
        #iterate through dictionary (keys are sample dirs in salmon dir)
        for key in counts_db: 
            #print("current key in dict: ", key)
            #change into salmon dir
            os.chdir(key)
            #print("you are now in: ", os.getcwd())
            
            #Reset add_to_line between each quant file
            global add_to_line
            add_to_line = ""
            
            #open quant file, add TPM and NumReads to list 
            with open("quant.sf", "r") as in_handle:
                #print("you are now in the quant.sf file") 
                for line in in_handle: 
                    line = line.rstrip()
                    sp_line = line.split("\t")
                    if sp_line[0].startswith(header): 
                        counts_db[key].append(sp_line[3]) #TPM
                        counts_db[key].append(sp_line[4]) #NumReads
                print("The header, sample, and list: ", header, key, counts_db[key])
       
            
            #Add quant info to the blastout info and write the one line of info to the output file
            #make addition of quant data to add to blastout line
            #global add_to_line 
            add_to_line = ",{0},{1},{2}".format(key, float(counts_db[key][0]), float(counts_db[key][1]))
            print("info from salmon to add to line: ", add_to_line)

            #make full line to write
            global write_line
            write_line = line_to_write + add_to_line
            print("line that will get written: ", write_line)
        
            #write csv line 
            with open("{0}/{1}.csv".format(cwd,args.d), "a", newline="") as out_handle:
                #split string into list to write to csv
                split_full_line = write_line.split(",")
                print("split line list: ", split_full_line)
                writer = csv.writer(out_handle, delimiter=",")
                writer.writerow(split_full_line)
                #out_handle.write(write_line)

            #Change out of sample dir to main salmon dir between each
            #print("after writing cwd: ", os.getcwd())
            os.chdir("..") 
            #print("after writing, changed to salmon dir: ", os.getcwd())
   
                  
    except IOError as err:
        print("problem reading or writing/appending file:", err)               
    
    
    #Change back to your original dir
    os.chdir(cwd)
    #print("changed back to original dir: ", os.getcwd())



#Part 2: Run Function: open blastout file and iterate through each line - on each line
#extract  header, prep the line to be written out and run the Function 
 
with open(args.a, "r") as in_handle:
    for line in in_handle:
        line = line.rstrip()
        sp_line = line.split("\t")
        #print(sp_line)
        
        header = sp_line[0].split("..")[0]
        #print(header)
        header_loc = sp_line[0].split("..")[1]
        
        line_to_write = "{0},{1},{2},{3},{4}".format(sp_line[0],sp_line[2],sp_line[1], header,header_loc)
        
        print("on header: ", header)
        print("using header in function for salmon lookup")
        
        #make global variable add_to_line
        add_to_line = ""
        
	    #Call second function using output of first function(the current header) and dir you specified 
        get_salmon_info(header, args.b)
 
 
print("Part 1 and 2 are complete!! Your output file is: ", args.d)
print("Moving to part 3: Rank and sort the output file")



#Part 3 Sort the output  
#open the csv file you made in Part 1 and 2 
reader = csv.reader(open("{0}.csv".format(args.d)), delimiter=",")

#skip the header
headers=next(reader)

#sort the data based on column 6 which is TPM 
sortedlist = sorted(reader, key=lambda row: float(row[6]), reverse=True)
print("The sorted list", sortedlist)   

#Write out the newly sorted file (header first then loop through the rest of the file)
with open("sorted_TPM_file.csv", "w", newline="") as out_handle:
    writer = csv.writer(out_handle, delimiter=",")
    writer.writerow(headers)

    for line in sortedlist:
        writer.writerow(line)

   
print("Program complete!! Your Ranked and Sorted output file is: sorted_TPM_file.csv")  


