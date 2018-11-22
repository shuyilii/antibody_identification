#!/usr/bin/python3

import pandas as pd
import sys

def get_protein_set(msgf_result):
    protein_list = []
    with open(msgf_result, 'r') as file:
        next(file)
        for line in file:
            element_list = line.rstrip().split('\t')
            protein = element_list[10]
            protein_list.append(protein)
        protein_set = set(protein_list)
    return protein_set

def get_protein_dict(msgf_result):
    protein_set = get_protein_set(msgf_result)
    protein_dict = {}
    for each in protein_set:
        protein_dict[each] = []
    msgf = pd.read_csv(msgf_result, sep = '\t')
    for i in range(len(msgf.index)):
        protein = msgf.iloc[i]['Protein']
        peptide = ''.join([i for i in msgf.iloc[i]['Peptide'][2:-2] if not i.isdigit()]).replace('+.','')
        protein_dict[protein].append(peptide)
    return protein_dict

def combine_dicts(dict1, dict2):
    final_dict = {}
    protein_set = set(set(dict1.keys()) & set(dict2.keys()))
    for each in protein_set:
        final_dict[each] = []
    for protein in dict1.keys():
        final_dict[protein].extend(dict1[protein])
    for protein in dict2.keys():
        final_dict[protein].extend(dict2[protein])
    return final_dict

combine_dict = get_protein_dict(sys.argv[1])
for i in range(2,len(sys.argv)):
    combine_dict = combine_dicts(combine_dict, get_protein_dict(sys.argv[i]))

#print(combine_dict)
