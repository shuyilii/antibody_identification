#!/usr/bin/python3
import sys

with open(sys.argv[1],'r') as fin:
    pro_dict = {}
    for line in fin:
        line = line.rstrip()
        if line.startswith('Protein'):
            protein = line.replace('Protein:', '')
        if line.startswith('Contigs'):
            contigs = line.replace('Contigs:', '').replace('[','').replace(']','').replace("'",'').split(', ')
            pro_dict[protein] = contigs
pep_list = []
for each in pro_dict.values():
    pep_list.extend(each)
pep_set = set(pep_list)
pep_dict = {}
for pep in pep_set:
    protein_list = []
    for protein in pro_dict:
        if pep in pro_dict[protein]:
            protein_list.append(protein)
    pep_dict[pep] = protein_list
print('Unique peptide number:'+ str(len(pep_dict)))
for each in pep_dict:
    print(each)
    print(str(pep_dict[each])+'\n')
