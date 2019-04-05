#!/usr/bin/python3

import sys
####import pandas as pd
from matplotlib import pyplot as plt

A_intens_tsv = sys.argv[1]
B_intens_tsv = sys.argv[2]

def intens(file):
    pep_dict = {}
    with open(file, 'r') as fin:
        next(fin)
        for line in fin:
            element = line.rstrip().split('\t')
            if element[1].find('IGH') > 0 and float(element[2]) < 0.18:
                pep = ''.join([i for i in element[0] if not i.isdigit()]).replace('+.','')
                intens = float(element[3])
                if pep not in pep_dict:
                    pep_dict[pep] = [intens]
                else:
                    pep_dict[pep].append(intens)
    max_pep_dict = {}
    for k,v in pep_dict.items():
        max_pep_dict[k] = max(v)
    return max_pep_dict

def rare_pep(file):
    rare_pep_list = []
    with open(file, 'r') as fin:
        next(fin)
        for line in fin:
            element = line.rstrip().split('\t')
            if element[1].find('IGH') > 0 and len(element[1].split(';')) < 4:
                pep = ''.join([i for i in element[0] if not i.isdigit()]).replace('+.','')
                rare_pep_list.append(pep)
    return rare_pep_list

A_intens_dict = intens(A_intens_tsv)
B_intens_dict = intens(B_intens_tsv)
A_rare_pep_list = rare_pep(A_intens_tsv)
B_rare_pep_list = rare_pep(B_intens_tsv)
rare_pep_list = set(A_rare_pep_list) & set(B_rare_pep_list)

intersect_dict = {}
for pep in A_intens_dict:
    if pep in B_intens_dict:
        intersect_dict[pep] = (A_intens_dict[pep],B_intens_dict[pep])

interest_pep = {}
d1 = [[],[]]
d2 = [[],[]]
d3 = [[],[]]
for k in intersect_dict:
    intens = intersect_dict[k]
    if intens[0] < 5e8:
        if k in rare_pep_list:
            fold = intens[0]/intens[1]
            if fold > 3:
                d1[0].append(intens[0])
                d1[1].append(intens[1])
                interest_pep[k] = fold
            else:
                d2[0].append(intens[0])
                d2[1].append(intens[1])
        else:
            d3[0].append(intens[0])
            d3[1].append(intens[1])

###visualization
data = (d3,d2,d1)
colors = ('k','g','r')
groups = ('not_interested','less_interested','most_interested')
for data, color, group in zip(data, colors, groups):
    x, y = data
    plt.scatter(x, y, c = color,label = group,s = 10)
plt.legend(loc = 'upper right')
for k,v in sorted(interest_pep.items(), key=lambda kv: kv[1], reverse = True):
    print(k + ': ' + str(v))
plt.xlabel('Php2', fontsize = 10)
plt.ylabel('IgG', fontsize = 10)
plt.title('scatter plot of peptide intensity')
plt.show()
