#!/usr/bin/python3
import sys

intens_file1 = sys.argv[1]
intens_file2 = sys.argv[2]

def parse(file):
    pep_dict = {}
    protein_list = []
    with open(file, 'r') as fin:
        next(fin)
        for line in fin:
            element = line.rstrip().split('\t')
            if element[1].find('IGH') > 0 and float(element[2]) < 0.18:
                protein_list.extend(element[1].split(';'))
                pep = ''.join([i for i in element[0] if not i.isdigit()]).replace('+.','')
                intens = float(element[3])
                if pep not in pep_dict:
                    pep_dict[pep] = [intens]
                else:
                    pep_dict[pep].append(intens)
    max_pep_dict = {}
    for k,v in pep_dict.items():
        max_pep_dict[k] = max(v)
    protein_list1 = [w[:-14] for w in protein_list]
    protein_set = set(protein_list1)
    return max_pep_dict,protein_set

def rank(intens_file):
    max_pep_dict = parse(intens_file)[0]
    protein_set = parse(intens_file)[1]
    assemble_dict = {}
    for each in protein_set:
        assemble_dict[each] = []
    with open(intens_file, 'r') as f1:
        next(f1)
        for line in f1:
            element = line.rstrip().split('\t')
            if element[1].find('IGH') > 0 and float(element[2]) < 0.18:
                protein = element[1].split(';')
                protein = [w[:-14] for w in protein]
                pep = ''.join([i for i in element[0] if not i.isdigit()]).replace('+.','')
                for each in protein:
                    assemble_dict[each].append(pep)
    intens_dict = {}
    for each in protein_set:
        intens_dict[each] = []
    for each in assemble_dict:
        assemble_dict[each] = set(assemble_dict[each])
        for pep in assemble_dict[each]:
            intens_dict[each].append(max_pep_dict[pep])
    for protein in intens_dict:
        intens_dict[protein] = sum(intens_dict[protein])
    r = {key: rank for rank, key in enumerate(sorted(set(intens_dict.values()), reverse=True), 1)}
    rank_dict = {k: r[v] for k,v in intens_dict.items()}
    return rank_dict
rank_dict1 = rank(intens_file1)
rank_dict2 = rank(intens_file2)
same_protein = set(rank_dict1.keys()) & set(rank_dict2.keys())
final_dict = {}
for protein in same_protein:
    final_dict[protein] = rank_dict1[protein]-rank_dict2[protein]
pos_dict = { k: v for k, v in final_dict.items() if v > 0 }
sorted_tuple = sorted(pos_dict.items(), key=lambda x: x[1], reverse=True)
for each in sorted_tuple:
    print(each[0] + '\t' + str(each[1]))
