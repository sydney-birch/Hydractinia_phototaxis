#! /usr/bin/env python3

#import modules
import argparse
import time

#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="input file: gene set file ")
parser.add_argument("-b", help="input file: gene_symbol_accid")

args = parser.parse_args()

t1 = time.time()

#key = seqid of interest (from blastout [1])
#value = [accid, OG, symbol]

tot_db = {}

#part 1 make db with seqid as key and a list with accid as the value
with open("missing_seqs.txt", "w") as out_handle:
    with open(args.a, "r") as in_handle:
        for line_a in in_handle:
            line_a = line_a.rstrip()
            sp_line_a = line_a.split(",")
            tot_db.setdefault(sp_line_a[1],[]) #set key to symbol
            tot_db[sp_line_a[1]].append(sp_line_a[0]) #adds entrez id

        print("part 1 dict: ", tot_db)
        print("number of keys: ", len(tot_db))

#part 2 add all human OGs to corresponding seqid of interest 
    with open(args.b, "r") as in_handle_2:
        for line_b in in_handle_2:
            line_b = line_b.rstrip()
            sp_line_b = line_b.split("\t")
            #print("split ogs: ", sp_line_b)
            if sp_line_b[0] in tot_db:
                if len(sp_line_b) == 1:
                    print("missing value(accid): ", sp_line_b[0])
                else: 
                    tot_db[sp_line_b[0]].append(sp_line_b[1])
            
        for key in tot_db:
            #print("len of list: ", len(tot_db[key]))
            if len(tot_db[key]) < 2: 
                #print("key: {0} and value: {1}".format(key, tot_db[key]))
                #print("entrez id: {0}".format(tot_db[key][0]))
                out_handle.write("{0}\t{1}\n".format(key, tot_db[key][0]))
