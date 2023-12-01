#! /usr/bin/env python3

#This script will clean up the new csv file you just created of a gene set from the Broad Institute GSEA (it gets rid of the trailing .. at the end of the discription and brackets).
#After this clean up step you can then move on to downloading the genes from the gene set in the next script.

#import modules
import argparse

#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input csv file with a list of gene symbols and descriptions - each on a single line")

args = parser.parse_args()

#Open the given input csv file and the mod file you are writing to
with open("{0}-mod".format(args.i), "w") as out_handle:
    with open(args.i, "r") as in_handle:
        #loop through each line of the file if there is a [ split on that or else split on ..  
        #and  write to the first poriton to a  new file 
        for line in in_handle: 
            line = line.rstrip()
            if "[" in line: 
                sp_line = line.split("[")
                #print("line split: ", sp_line)
                write_out = sp_line[0] + "\n"
                out_handle.write(write_out)
            elif ".." in line: 
                sp_line = line.split("..")
                #print("split line: ", sp_line)
                write_out = sp_line[0] + "\n"
                out_handle.write(write_out)
        print("Your Modified File has been Created")
