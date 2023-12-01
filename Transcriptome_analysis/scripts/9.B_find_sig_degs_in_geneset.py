#! /usr/bin/env python3

#The goal of this script is to identify significant Actinula Transcripts in the given gene set. This script will compare the gene set output file from R (made at step 8)
#to the edgeR output (you will need to run this for how ever many edgeR comparison files you have). It creates a dictionary based on the headers from the step 8 file
#and will see if any of those headers are in the edgeR output file - if there is a significant gene present it will record all of the info from the edgeR file and will write
#it out to an output file. 

#import modules
import argparse
import time

#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="input file: gene set from R (file with gene_acc OG symbol Actinula_seqid)")
parser.add_argument("-b", help="input file: DEG list - made from edgeR")
parser.add_argument("-c", help="name of output file (the stage comparison)")

args = parser.parse_args()

#I want to go through the input file of our genes of interest, add those as keys to a dictionary
#then iterate through the DEG list and collect data for any headers that match and report that back 

t1 = time.time()

#dict set up: 
#key = header (pos 4) in file a (in file b pos 0)
#value = [OG,symbol, seqid, logFC, pvalue]
	# 0    1       2     3     4 
	
tot_db = {}
#part 1 populate dictionary with headers of interest from gene set, and open output file and prep it
with open("{0}_DEGs.txt".format(args.c), "w") as out_handle:
    out_handle.write("Header\tOG\tSymbol\tlogFC\tlogCPM\tPvalue\n")
    with open(args.a, "r") as in_handle:
        for line_a in in_handle:
            line_a = line_a.rstrip()
            sp_line_a = line_a.split("\t")
            #print(sp_line_a)
            #sp_line_a[3].strip()
            tot_db.setdefault(sp_line_a[4].strip(),[]) #set key to header and make list as value
            tot_db[sp_line_a[4]].append(sp_line_a[2].strip()) #add OG
            tot_db[sp_line_a[4]].append(sp_line_a[3].strip()) #add symbol
            
        print("db after 1st file: ", tot_db, len(tot_db))   

        #if there are dups in the list remove here (in future make this better)
        for key in tot_db: 
            print("length of cur val", len(tot_db[key]))
            if len(tot_db[key]) == 16:
                tot_db[key].pop()
                #print("length of list after removing last item: ", len(tot_db[key]))
                tot_db[key].pop()
                tot_db[key].pop()
                tot_db[key].pop()
                tot_db[key].pop()
                tot_db[key].pop()
                tot_db[key].pop()
                tot_db[key].pop()
                tot_db[key].pop()
                tot_db[key].pop()
                tot_db[key].pop()
                tot_db[key].pop()
                tot_db[key].pop()
                tot_db[key].pop()
                #print("length of list after removing last item second time: ", len(tot_db[key]))
            elif len(tot_db[key]) == 8:
                tot_db[key].pop()
                #print("length of list after removing last item: ", len(tot_db[key]))
                tot_db[key].pop()
                tot_db[key].pop()
       	       	tot_db[key].pop()
       	       	tot_db[key].pop()
       	       	tot_db[key].pop()
                #print("length of list after removing last item second time: ", len(tot_db[key]))
            elif len(tot_db[key]) == 4: 
                tot_db[key].pop()
                #print("length of list after removing last item: ", len(tot_db[key]))
                tot_db[key].pop()
                #print("length of list after removing last item second time: ", len(tot_db[key]))
        print("db after editing: ", tot_db)
        print("length of db: ", len(tot_db))
 
        print("now checking for matches in edgeR file ")
            
#part 2 add logFC, logCPM, and P value to list in dict if header matches  
    with open(args.b, "r") as in_handle_2:
        for line_b in in_handle_2:
            line_b = line_b.rstrip()
            sp_line_b = line_b.split("\t")  
            sp_line_b[0].strip()
            #print("header in r file", sp_line_b[0])
            #print("logfc: ", sp_line_b[1])

            if sp_line_b[0].strip() in tot_db: 
                tot_db[sp_line_b[0]].append(sp_line_b[1]) #add logFC
                tot_db[sp_line_b[0]].append(sp_line_b[2]) #add logCPM
                tot_db[sp_line_b[0]].append(sp_line_b[3]) #add pval
                print("match: {0} logFC: {1}".format(sp_line_b[0], sp_line_b[1]))
         
        print("fully populated db: ", tot_db, len(tot_db))
         
#part 3: write to file
    for key in tot_db: 
        print("Writing to file, cur value length: ", len(tot_db[key]))
        #print("len with 0", len(tot_db[key][0]))
        #print(tot_db[key])
        if len(tot_db[key]) == 5: 
            print("key: ", key)
            print("logFC: ", tot_db[key][2])
            out_handle.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(key, tot_db[key][0], tot_db[key][1], tot_db[key][2], tot_db[key][3], tot_db[key][4]))                           
             
#time program
t2 = time.time()
secs_time = t2-t1
int(secs_time)
tot_time = secs_time/60

print("Total running time: {0} minutes ".format(tot_time)) #in minutes
