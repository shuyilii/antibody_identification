#!/usr/bin/python3

import pandas as pd
from itertools import combinations
import sys
import matplotlib.pyplot as plt

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

def get_ref_dict(ref_database):
    with open(ref_database, 'r') as database:
        ref_database_dict = {}
        for line in database:
            line = line.rstrip()
            if line.startswith('>'):
                protein_ID = line.replace('>','')
            else:
                ref_database_dict[protein_ID] = line
    return ref_database_dict

def get_mapping_dict(combine_dict,ref_database):
    ref_database_dict = get_ref_dict(ref_database)
    for each in combine_dict:
        combine_dict[each] = set(combine_dict[each])
    mapping_dict = {}
    for protein in ref_database_dict:
        temp_list = []
        for each in combine_dict[protein]:
            mapping_pos = (ref_database_dict[protein].find(each) + 1, ref_database_dict[protein].find(each) + len(each))
            temp_list.append(mapping_pos)
        mapping_dict[protein] = temp_list
    return mapping_dict

def merge_reads(position_list):
    saved = list(position_list[0])
    for st, en in sorted([sorted(t) for t in position_list]):
        if st <= saved[1] + 1:
            saved[1] = max(saved[1], en)
        else:
            yield tuple(saved)
            saved[0] = st
            saved[1] = en
    yield tuple(saved)

def visualization(each_mapping_list, ref_database_dict, each_merged_list, protein):
    max_length = len(max(list(ref_database_dict.values()),key = len))
    reads_num = len(each_mapping_list) + len(each_merged_list)
    ax = plt.axes(xlim = (0,max_length + 50),ylim = (0,reads_num + 2),
    xticks = range(0, max_length + 50, 40),title = protein)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.arrow(1,1,len(ref_database_dict[protein])-1,0, head_width=0.05, head_length=0.1, fc='r', ec='r')
    i = 2
    for merged_pos in each_merged_list:
        ax.arrow(merged_pos[0], i, merged_pos[1]-merged_pos[0], 0, head_width=0.05, head_length=0.1, fc='g', ec='g')
        i += 1
    for mapping_pos in each_mapping_list:
        ax.arrow(mapping_pos[0], i, mapping_pos[1]-mapping_pos[0], 0, head_width=0.05, head_length=0.1, fc='k', ec='k')
        i += 1
    plt.show()

combine_dict = get_protein_dict(sys.argv[2])
for i in range(3,len(sys.argv)):
    combine_dict = combine_dicts(combine_dict, get_protein_dict(sys.argv[i]))

mapping_dict = get_mapping_dict(combine_dict,sys.argv[1])

merged_dict = {}
for protein in mapping_dict:
    sorted_list = sorted(mapping_dict[protein], key=lambda tup: tup[0])
    merged_dict[protein] = list(merge_reads(sorted_list))

ref_database_dict = get_ref_dict(sys.argv[1])
for protein in ref_database_dict:
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
    print('Protein:' + protein)
    print('Coverage:' + str(coverage) + '%')
    print('Contigs:' + str(temp_seq_list))
    print('Position:' + str(pos_list) + '\n')

protein = input("Which protein do you want to visualize?:")
visualization(mapping_dict[protein],ref_database_dict, merged_dict[protein], protein)
