#! /usr/bin/env python3

#This script will take in your handmade summary file of sig DEGs (file that you saved in a text file from 9.B_find_sig_degs_in_geneset.py that has each comparisons output). 
#It will take all of the sig DEG headers and put them in the format of: "header_1","header_2","header_3"..... so you can import this string into R to make heatmaps. 
#copy the output in the output file and add it to your summary file.  

#import modules
import argparse
import time

#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="input file: summary file of found sig DEGs")
parser.add_argument("-b", help="name of output file - name of gene set")

args = parser.parse_args()

t1 = time.time()

	
tot_list = []
#part 1 populate list with headers and open output file to write to 
with open("{0}-list_of_sig_DEGs.txt".format(args.b), "w") as out_handle:
    with open(args.a, "r") as in_handle:
        for line in in_handle:
            line = line.rstrip()
            sp_line = line.split(" ")
            #print(sp_line)
            
            #add headers to the list which are in the 0 position
            if "_t." in sp_line[0]:
                tot_list.append(sp_line[0])
            
        print("Full list populated: ", tot_list, len(tot_list))   

#part 2 join list in correct format for R and write to output file 
new_string = "\",\"".join(tot_list)
print(new_string)

with open("{0}-list_of_sig_DEGs.txt".format(args.b), "w") as out_handle:
    out_handle.write("\"{0}\"".format(new_string))


#time program
t2 = time.time()
secs_time = t2-t1
int(secs_time)
tot_time = secs_time/60

print("Total running time: {0} minutes ".format(tot_time)) #in minutes
