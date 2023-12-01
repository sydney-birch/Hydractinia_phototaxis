#! /usr/bin/env python3

#Dont forget to module load linuxbrew/colsa 

#This script renames the headers after running transdecoder. 
#A typical header looks like: >Gene.1::Ec_actinula_t.1::g.1::m.1 type:complete len:236 gc:universal Ec_actinula_t.1:2328-1621(-)
#This script will use the last segment of the transdecoder header but without the direction (-) or (+) and without the colon
#so the resulting header will be:  >Ec_actinula_t.1..2328-1621. so we can tell which original transcript it came from and what section of that transcript
#The output file will have the -mod.fa tag at the end of the input fasta file 

#import modules
import argparse

#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="input transdecoder fasta file")

args = parser.parse_args()
#open the input and the output file so you can read and write simultaneously
try: 
    with open(args.a, "r") as in_handel: 
        with open("{0}-mod.fa".format(args.a), "w") as out_handel: 
            #go through each line of the file if the file is a header, split the header and rejoin it to get the new header 
            for line in in_handel: 
                if line.startswith(">"): 
                    sp_line = line.split(" ")
                    #print("first split: ", sp_line)
                    sp_line_2=sp_line[4].split("(")
                    #print("second split: ", sp_line_2)
                    sp_line_3= sp_line_2[0].split(":")
                    join_sp= "..".join(sp_line_3)
                    #print("new line after join", join_sp)                

                    new_line = ">{0}".format(join_sp)+ "\n"
                    #print("new line: ", new_line)
                    out_handel.write(new_line)
                    
                else: 
                    out_handel.write(line)

except IOError as err:
    print("problem reading or writing/appending file:", err)
