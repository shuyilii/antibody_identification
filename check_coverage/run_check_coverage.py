#!/usr/bin/python3

import argparse
from check_coverage import *

file_num = len(sys.argv) - 10
parser = argparse.ArgumentParser(description="sample usage: run_check_coverage.py -r ref_database.fasta -i data1.tsv data2.tsv data3.tsv ... -c 0.01")
parser.add_argument("-r", "--ref", help = "Reference fasta file")
parser.add_argument("-i", "--input", help = "Input msgf+ tsv file", nargs=file_num)
parser.add_argument("-c1", "--cutoff1", help = "The strict cutoff threshold of Q-value",type = float)
parser.add_argument("-c2", "--cutoff2", help = "The loose cutoff threshold of Q-value",type = float)
parser.add_argument("-fr_pos", "--fr_pos", help = "The position information of framework and cdr region")
args = parser.parse_args()

if len(sys.argv) < 11:
    print('usage: run_check_coverage.py [-h] [-r REF] [-i] [-c1 CUTOFF] [-c2 CUTOFF2] [-fr_pos FR_POS]')
    print('use "run_check_coverage.py -h" for more information')
    sys.exit()

combine_dict = get_protein_dict(args.input[0], args.cutoff1, args.cutoff2)
for i in range(1,len(args.input)):
    combine_dict = combine_dicts(combine_dict, get_protein_dict(args.input[i], args.cutoff1, args.cutoff2))

mapping_dict = get_mapping_dict(combine_dict,args.ref)

merged_dict = {}
for protein in mapping_dict:
    sorted_list = sorted(mapping_dict[protein], key=lambda tup: tup[0])
    merged_dict[protein] = list(merge_reads(sorted_list))

ref_database_dict = get_ref_dict(args.ref)
coverage_dict = {}
for protein in ref_database_dict:
    if protein in merged_dict:
        full_length = len(ref_database_dict[protein])
        pos_list = merged_dict[protein]
        temp_pos_list = []
        temp_seq_list = []
        for each in pos_list:
            temp_length = each[1] - each[0] + 1
            temp_pos_list.append(temp_length)
            temp_seq_list.append(ref_database_dict[protein][each[0]-1:each[1]])
        mapping_length = sum(temp_pos_list)
        coverage = (mapping_length/full_length)*100
        coverage_dict[protein] = coverage
        if len(temp_seq_list) > 2:
            print('Protein:' + protein)
            print('Coverage:' + str(coverage) + '%')
            print('Contigs:' + str(temp_seq_list))
            print('Position:' + str(pos_list) + '\n')

while 1:
    protein = input("Which protein do you want to visualize?(print exit if you want to exit):")
    if protein == 'exit':
        sys.exit()
    else:
        file_dict = {}
        for file in args.input:
            if protein in get_protein_dict(file, args.cutoff1, args.cutoff2):
                temp_list = []
                file_dict[file] = get_protein_dict(file, args.cutoff1, args.cutoff2)[protein]
            else:
                file_dict[file] = []

        file_mapping_dict = {}
        for each in file_dict:
            temp_list = []
            for pep in file_dict[each]:
                mapping_pos = (ref_database_dict[protein].find(pep) + 1, ref_database_dict[protein].find(pep) + len(pep))
                temp_list.append(mapping_pos)
            file_mapping_dict[each] = set(temp_list)

        pos_dict = {}
        with open(args.fr_pos, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    ID = line.rstrip().replace('>','')
                else:
                    pos_dict[ID] = line.rstrip().split(' ')

        visualization(file_mapping_dict, ref_database_dict, pos_dict,merged_dict[protein], protein, coverage_dict[protein])
