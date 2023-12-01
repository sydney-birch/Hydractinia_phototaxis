#! /usr/bin/env python3

#Don't forget to module load linuxbrew/colsa

#Use if you have 150-500 genes in gene set (Part 1) This script will use the csv file you just created and cleaned up to search for and download the verified protien 
#sequences of all of the genes in your selected gene set. You will need to change your email address for the NCBI search on line 37. 
#This script will do steps: 1) initial search to get accids, and 2) make a dict of gene symbols and accids. 
#The next step will be completed in the next script which is step 3) use dict to search and make fasta.  


#import modules
import argparse
import Bio.Entrez
import time
import json

#create an instance of Argument Parser and add positional argument
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input csv file with a list of entrez ids, gene symbols and discriptions - each gene on a single line")
#parser.add_argument("-o", help="output fasta file name")

args = parser.parse_args()

###This is get_entrez_fastas-v5 split up for larger data sets - part 1

#Part 1 - do intial search to get the accids
with open("gene_tables.txt", "w") as out_handle:
    with open("summary_info.txt", "w") as out_handle_2:
        with open("tmp_dict_for_searches.txt", "w") as out_handle_3:
            with open(args.i, "r") as in_handle:
                #Loop through input file of gene symbols - for each line(gene) search ID of interest in Homo sapiens
                for line in in_handle:
                    line = line.rstrip()
                    sp_line = line.split(",")  #entrez id   symbol   description
                    print("current search term: ", sp_line)
                    time.sleep(3)
                    Bio.Entrez.email="sjb1061@wildcats.unh.edu"
                    #handle = Bio.Entrez.esearch(db="protein", term="({0}[text word]) AND Homo sapiens[organism] AND {1}[title]".format(sp_line[0], sp_line[1]))
                    handle = Bio.Entrez.esearch(db="gene", term="({0}[uid]) AND Homo sapiens[organism]".format(sp_line[0])) #taking entrez id
                    results = Bio.Entrez.read(handle)
                    handle.close()
                    #print("not in if statement: ", results)

                    if int(results["Count"]) > 0:
                        #print("results dict: ", results) #debug
                        #put to sleep to help with online search use a try loop
                        time.sleep(4)
                        try:
                            #get data from search - take the first hit from the search
                            handle = Bio.Entrez.efetch(db="gene", id=results["IdList"][0], rettype="gene_table", retmode="text")
                            handle2 = Bio.Entrez.efetch(db="gene", id=results["IdList"][0], rettype="acc", retmode="text")
                            fasta_data=handle.read() #gene table
                            fasta_data2=handle2.read() #acc - writing summary info
                            handle.close()

                            #write to output file
                            #print("fasta data: ", fasta_data) #debug
                            #print("fasta data 2: ", fasta_data2)
                            fasta_data2 = fasta_data2.rstrip()
                            out_handle_2.write("{0}\t{1}\n".format(sp_line[1],fasta_data2))
                            out_handle.write(fasta_data)

                        except HTTPError:
                            print("Need more time between search submissions")

                    else:
                        print("search for this by hand: ", line)
                        write_out = line + "\n"
                        out_handle_3.write(write_out)

print("part 1 complete - moving to part 2")

#part 2 - populate a dictionary with the accids as the values and the symbols as the keys
acc_db = {}

with open("gene_tables.txt", "r") as in_handle:
    cur_key = ""
    print("cur key at start", cur_key)
    for line in in_handle:
        if "[Homo sapiens]" in line:
            sp_header = line.split(" ")
            acc_db.setdefault(sp_header[0], [])
            cur_key = sp_header[0]
            #print("current key", cur_key)

        if line.startswith("protein"):
            sp_line = line.split(" ")
            #print("split line: ", sp_line)
            for item in sp_line:
                if "P_" in item:
                    acc_db[cur_key].append(item)

            #if sp_line[1].startswith("NP_"):
                #print("acc id starts with NP: ", sp_line[1])
                #print("adding to cur key: ", cur_key)
                #acc_db[cur_key] = sp_line[1]
            #else:
                #print("acc id on second spot: ", sp_line[2])
                #print("adding to cur key: ", cur_key)
                #acc_db[cur_key] = sp_line[2]

    print("tot db and length: ", acc_db, len(acc_db))
    with open("tmp_dict_for_searches.txt", "w") as out_handle_3:
        out_handle_3.write(json.dumps(acc_db))
    print("dictionary has been written to new file: tmp_dict_for_searches.txt move to next script")
    
