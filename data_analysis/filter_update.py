#!/usr/bin/python3

import sys
from colors import *
####input files
protein_result = sys.argv[1]
uniq_pep = sys.argv[2]
####input args
fil_num = int(input('The rare peptide should be in no more than how many antibodies:'))
cov_thres = float(input('The coverage threshold:'))
cdr3 = input('If only pick those antibodies that have peptides in cdr3[Y or N]:')

with open(uniq_pep, 'r') as f1:
    next(f1)
    pep_dict = {}
    for line in f1:
        line = line.rstrip()
        if line.startswith('peptide'):
            pep = line.replace('peptide:','')
        if line.startswith('protein'):
            pro = line.replace('protein:','').replace('[','').replace(']','').replace("'",'').split(', ')
            pep_dict[pep] = pro

filtered_pro_id_dict = dict((k, v) for k, v in pep_dict.items() if len(v) <= fil_num)

def colored_seq(ref_seq,rare_pep,contigs_hit_pos):
    rare_pep_pos = []
    for pep in rare_pep:
        start_pos = ref_seq.find(pep)
        end_pos = start_pos + len(pep)
        rare_pep_pos.append((start_pos, end_pos))
    common_pep_pos = [x for x in contigs_hit_pos if x not in rare_pep_pos]
    colored_seq = ''
    for i in range(len(ref_seq)):
        if any(start <= i < end for start,end in rare_pep_pos):
            colored_seq += red(ref_seq[i])
        elif any(start <= i < end for start,end in common_pep_pos):
            colored_seq += green(ref_seq[i])
        else:
            colored_seq += ref_seq[i]
    return colored_seq

with open(protein_result, 'r') as f2:
    info_dict = {}
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
            if line.startswith('Position'):
                pos = line.replace('Position:','').replace('[','').replace(']','').replace('(','').replace(')','').split(',')
                temp_pos_list = list(zip(*[iter(pos)] * 2))
                contigs_hit_pos = []
                for pos in temp_pos_list:
                    contigs_hit_pos.append((int(pos[0])-1,int(pos[1])))
            elif line.startswith('Coverage'):
                coverage = float(line.replace('Coverage:','').replace('%',''))/100
                info_dict[protein] = [coverage]
            elif line.startswith('Ref_seq'):
                ref_seq = line.replace('Ref_seq:','')
                info_dict[protein].append([ref_seq,rare_pep,contigs_hit_pos])
            elif line.startswith('Region'):
                info_dict[protein].append(line)

for protein in info_dict:
    if info_dict[protein][0] > cov_thres:
        pro_style = blue(protein, style='bold+underline')
        print('Protein:'+ pro_style)
        ref_seq = info_dict[protein][1][0]
        rare_pep = info_dict[protein][1][1]
        contigs_hit_pos = info_dict[protein][1][2]
        print('Ref_seq:' + colored_seq(ref_seq,rare_pep,contigs_hit_pos))
        print(contigs)
        print(info_dict[protein][2]+'\n')
