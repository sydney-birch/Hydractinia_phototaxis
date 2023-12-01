#! /usr/bin/env python3

#Dont forget to module load linuxbrew/colsa 

#This script takes in 2 arguments, the fasta file that you want to change the headers in and what you want the header to change to and it writes 
#a new modified fasta file with the altered headers. The format of the new headers will be: >argument_b_t.# 
#where the # is the position of the header (the scripts counts the headers) and the t represents transcript. 
#For the actinula transcriptome it will look like Ec_actinula_t.1 : ./3.A_rename_fa_headers.py -a actinula_total.ORP.fa -b Ec_actinula 

#import modules
import argparse

#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="input fasta file")
parser.add_argument("-b", help="header name")

args = parser.parse_args()
#open the fasta file and the new file you are making which will have -mod.fa at the end of the input fasta file 
try: 
    with open(args.a, "r") as in_handel: 
        with open("{0}-mod.fa".format(args.a), "w") as out_handel: 
            
            #count the number of headers - if a header shows up write new header out, if the line is not a header write that line out to the file
            count = 1
            for line in in_handel: 
                if line.startswith(">"): 
                    new_line = ">{0}_t.{1}".format(args.b, count)+"\n"
                    #print("new line: ", new_line)
                    out_handel.write(new_line)
                    count +=1
                else: 
                    out_handel.write(line)

except IOError as err:
    print("problem reading or writing/appending file:", err)
