#! /usr/bin/env python3

#Don't forget to module load linuxbrew/colsa

#Use if you have 1-150 genes in gene set. This script will use the csv file you just created and cleaned up to search for and download the verified protien 
#sequences of all of the genes in your selected gene set. It will output 2 files that you will use in later step a FASTA file and a gene_symbol_accid file which will be 2 columns (tab delimited)
#the first column with have the gene symbol and the 2nd will have the corresponding NCBI verified Accession ID. You will also get 3 other ouput files: gene_tables, 
#summary_info, and search_by_hand. You will need to change your email address for the NCBI search on line 33 and 128. This script has 3 steps: 1) initial search to get accids, 2) make a dict of gene symbols and accids, 3) use dict to search and make fasta

#import modules
import argparse
import Bio.Entrez
import time

#create an instance of Argument Parser and add positional argument
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input cleaned up csv file with a list of gene symbols and info - each on a single line")
parser.add_argument("-o", help="output fasta file name")

args = parser.parse_args()
t1 = time.time()
#Part 1 - do intial search to get the accids
with open("gene_tables.txt", "w") as out_handle:
    with open("summary_info.txt", "w") as out_handle_2:
        with open("search_by_hand.txt", "w") as out_handle_3:
            with open(args.i, "r") as in_handle:
                #Loop through input file of entrez ids, gene symbols, and descriptions - for each line (gene) search ID of interest in Homo sapiens
                for line in in_handle:
                    line = line.rstrip()
                    sp_line = line.split(",")  #entrez id   symbol   description
                    print("current search term: ", sp_line)
                    time.sleep(4)
                    Bio.Entrez.email="sjb1061@wildcats.unh.edu"
                    #handle = Bio.Entrez.esearch(db="protein", term="({0}[text word]) AND Homo sapiens[organism] AND {1}[title]".format(sp_line[0], sp_line[1]))
                    handle = Bio.Entrez.esearch(db="gene", term="({0}[uid]) AND Homo sapiens[organism]".format(sp_line[0])) #taking entrez id
                    results = Bio.Entrez.read(handle)
                    handle.close()
                    #print("not in if statement: ", results)

                    if int(results["Count"]) > 0:
                        #print("results dict: ", results) #debug
                        #put to sleep to help with online search use a try loop
                        time.sleep(5)
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
time.sleep(10)

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

time.sleep(10)
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
                            #print("value is (should be NP with no comma): ", val)
                        else:
                            val = accid
                            #print("value is (should be NP): ", val)
            else:
                if "," in acc_db[key][0]:
                    sp_acc = acc_db[key][0].split(",")
                    val = sp_acc[0]
                    #print("value is (should be NP with no comma): ", val)
                else:
                    val = acc_db[key][0]
                    #print("value: ", val)
            print("current search term: ", acc_db[key][0])
            time.sleep(5)

            Bio.Entrez.email="sjb1061@wildcats.unh.edu"
            handle = Bio.Entrez.esearch(db="protein", term="({0}[accession])".format(val)) #search the current accession id
            results = Bio.Entrez.read(handle)
            handle.close()

            if int(results["Count"]) > 0:
                #print("results dict: ", results) #debug
                #put to sleep to help with online search use a try loop
                time.sleep(5)
                try:
                    #get data from search - take the first hit from the search
                    handle = Bio.Entrez.efetch(db="protein", id=results["IdList"][0], rettype="fasta", retmode="text")
                    fasta_data=handle.read() #fasta for current accid
                    handle.close()

                    #write to 2 output files
                    #print("fasta data: ", fasta_data) #debug should be fasta
                    #print("acc id: ", val)
                    #print("gene symbol: ", key)
                    out_handle.write(fasta_data)
                    out_handle_2.write("{0}\t{1}\n".format(key, val))

                except HTTPError:
                    print("Need more time between search submissions")

        print("your fasta file has been written!")

#time program
t2 = time.time()
secs_time = t2-t1
int(secs_time)
tot_time = secs_time/60

print("Total running time: {0} minutes ".format(tot_time)) #in minutes
