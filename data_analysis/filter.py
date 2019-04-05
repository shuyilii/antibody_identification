#!/usr/bin/python3
import sys
from itertools import islice

protein_result = sys.argv[1]
uniq_pep = sys.argv[2]
fil_num = sys.argv[3]

with open(uniq_pep, 'r') as f1:
    next(f1)
    pep_dict = {}
    for line in f1:
        line = line.rstrip()
        if line.startswith('peptide'):
            pep = line.rstrip().replace('peptide:','')
        if line.startswith('protein'):
            pro = line.rstrip().replace('protein:','').replace('[','').replace(']','').replace("'",'').split(', ')
            pep_dict[pep] = pro
filtered_pro_id_dict = dict((k, v) for k, v in pep_dict.items() if len(v) <= int(fil_num))
with open(protein_result, 'r') as f2:
    for line in f2:
        line = line.rstrip()
        if line.startswith('Protein'):
            hit = False
            pro = line.replace('Protein:','')
            rare_pep = []
            for each in filtered_pro_id_dict:
                if pro in filtered_pro_id_dict[each]:
                    rare_pep.append(each)
                    protein = pro
                    hit = True
            if hit:
                print('Protein:' + protein)
                print('Rare peptide:' + str(rare_pep))
                print(''.join(islice(f2, 6)))
