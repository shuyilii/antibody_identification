#!/usr/bin/python3

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

def get_protein_dict(msgf_result,cutoff1,cutoff2):
    protein_set = get_protein_set(msgf_result)
    protein_dict = {}
    protein_list_passed_strict_cutoff = []
    for each in protein_set:
        protein_dict[each] = []
    ####too slow for using
    #msgf = pd.read_csv(msgf_result, sep = '\t')
    # for i in range(len(msgf.index)):
    #     protein = msgf.iloc[i]['Protein']
    #     Q_value = msgf.iloc[i]['QValue']
    #     peptide = ''.join([i for i in msgf.iloc[i]['Peptide'][2:-2] if not i.isdigit()]).replace('+.','')
    #     if Q_value <= cutoff:
    #         protein_dict[protein].append(peptide)
    with open(msgf_result, 'r') as msgf:
        next(msgf)
        for line in msgf:
            line = line.rstrip()
            element_list = line.split('\t')
            protein = element_list[10]
            Q_value = float(element_list[15])
            if Q_value <= cutoff1:
                protein_list_passed_strict_cutoff.append(protein)
    with open(msgf_result, 'r') as msgf:
        next(msgf)
        for line in msgf:
            line = line.rstrip()
            element_list = line.split('\t')
            peptide = ''.join([i for i in element_list[9][2:-2] if not i.isdigit()]).replace('+.','')
            protein = element_list[10]
            Q_value = float(element_list[15])
            if protein in protein_list_passed_strict_cutoff:
                if Q_value <= cutoff2:
                    protein_dict[protein].append(peptide)
    final_protein_dict = {k:v for k,v in protein_dict.items() if v}
    return final_protein_dict

def combine_dicts(dict1, dict2):
    final_dict = {}
    protein_set = set(set(dict1.keys()) | set(dict2.keys()))
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
        if protein in combine_dict.keys():
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

def visualization(file_mapping_dict, ref_database_dict, pos_dict, each_merged_list, protein, coverage):
    max_length = len(ref_database_dict[protein])
    file_num = len(file_mapping_dict)
    color_list = ['burlywood','c','m','r','b','y','chartreuse','burlywood','c','m','r','b','y','chartreuse']
    ax = plt.axes(xlim = (-10 , max_length + 5),ylim = (-file_num - 2,file_num + 2 ),title = protein)
    ax.axis('off')
    ax.arrow(1+10,1,max_length-1 , 0 , head_width=0, head_length=0, fc='k', ec='k',width=0.2)
    plt.text(-40, 1, 'ref_seq' , fontsize=10)
    ax.arrow(1+10,2,max_length-1 , 0 , head_width=0, head_length=0, fc='k', ec='k',width=0.1)
    n = 0
    for pos in pos_dict[protein][0:6]:
        color = color_list[n]
        ax.arrow(int(pos)+10, 2, int(pos_dict[protein][6])-int(pos), 0, head_width=0, head_length=0, fc=color, ec=color,width=0.1)
        n += 1
    for merged_pos in each_merged_list:
        ax.arrow(merged_pos[0]+10, 1, merged_pos[1] - merged_pos[0], 0, head_width=0, head_length=0, fc='g', ec='g',width=0.2)
    i = 3
    temp_dict = {}
    for file in file_mapping_dict:
        color = color_list[i-2]
        if len(file_mapping_dict[file]) != 0:
            for mapping_pos in file_mapping_dict[file]:
                ax.arrow(mapping_pos[0]+10, i, mapping_pos[1] - mapping_pos[0], 0, head_width=0, head_length=0, fc=color, ec=color,width=0.1)
                plt.text(-40, i, str(i-2) , fontsize=10)
        else:
            plt.text(-40, i, str(i-2) , fontsize=10)
        temp_dict[str(i-2)] = file
        i +=1
    text_coverage = 'Coverage = ' + str(round(coverage,2)) + '%'
    a = file_num - 2
    for each in temp_dict:
        file_name = 'Dataset ' + str(each) + ':' + temp_dict[each]
        plt.text(-10, -file_num+a, file_name , fontsize=10)
        a += -1
    plt.text(max_length-80, -0.5, text_coverage , fontsize=10)
    plt.show()
