#! /usr/bin/env python3

#Don't forget to module load linuxbrew/colsa

#Use if you have 1-150 genes in gene set (Part 2) This script continues from part 1 - it downloads the verified protien 
#sequences of all of the genes in your selected gene set. You will need to change your email address for the NCBI search on line 64. 
#The previous script completed steps: 1) initial search to get accids, and 2) make a dict of gene symbols and accids. 
#This script runs step 3) use dict to search and make fasta. It will use the temporary file created from the first script and will generate the FASTA file. 

#import modules
import argparse
import Bio.Entrez
import time
import json

#create an instance of Argument Parser and add positional argument
parser = argparse.ArgumentParser()
#parser.add_argument("-i", help="input file with a list of gene symbols - each on a single line")
parser.add_argument("-o", help="output fasta file name")

args = parser.parse_args()

###This is get_entrez_fastas-v5 split up for larger data sets - part 2

#open dictionary from temp file and make it a global variable

with open("tmp_dict_for_searches.txt", "r") as in_handle_2:
    acc_db = json.loads(in_handle_2.read())

print("length of acc_db after opening file: ", len(acc_db))
print("acc_db: ", acc_db)

#Part 3: get the fastas for the accids
with open(args.o, "w") as out_handle:
    with open("gene_symbol_accid", "w") as out_handle_2:
        for key in acc_db:
            val = ""
            if len(acc_db[key]) > 1:
                for accid in acc_db[key]:
                    if accid.startswith("NP_"):
                        if "," in accid:
                            sp_acc = accid.split(",")
                            val = sp_acc[0]
                            print("value is (should be NP with no comma): ", val)
                        else:
                            val = accid
                            print("value is (should be NP): ", val)
            else:
                print("length of value in dict(should be 1): ", len(acc_db[key]))
                print("the list: ", acc_db[key])
                if len(acc_db[key]) == 0: 
                    with open ("search_by_hand.txt", "w") as out_handle_3:
                        out_handle_3.write("{0}\n".format(key))
                elif "," in acc_db[key][0]:
                    sp_acc = acc_db[key][0].split(",")
                    val = sp_acc[0]
                    print("value is (should be NP with no comma): ", val)
                else:
                    val = acc_db[key][0]
                    print("value: ", val)
            #print("current search term: ", acc_db[key][0])
            time.sleep(3)

            Bio.Entrez.email="sjb1061@wildcats.unh.edu"
            handle = Bio.Entrez.esearch(db="protein", term="({0}[accession])".format(val)) #search the current accession id
            results = Bio.Entrez.read(handle)
            handle.close()

            if int(results["Count"]) > 0:
                #print("results dict: ", results) #debug
                #put to sleep to help with online search use a try loop
                time.sleep(4)
                try:
                    #get data from search - take the first hit from the search
                    handle = Bio.Entrez.efetch(db="protein", id=results["IdList"][0], rettype="fasta", retmode="text")
                    fasta_data=handle.read() #fasta for current accid
                    handle.close()

                    #write to 2 output files
                    #print("fasta data: ", fasta_data) #debug should be fasta
                    print("acc id: ", val)
                    print("gene symbol: ", key)
                    out_handle.write(fasta_data)
                    out_handle_2.write("{0}\t{1}\n".format(key, val))

                except HTTPError:
                    print("Need more time between search submissions")

        print("your fasta file has been written!")
