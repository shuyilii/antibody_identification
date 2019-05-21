#!/usr/bin/python3

import sys
from colors import *
####input files
protein_result = sys.argv[1]
uniq_pep = sys.argv[2]
pos_file = sys.argv[3]
####input args
fil_num = int(input('The rare peptide should be in no more than how many antibodies:'))
cov_thres = float(input('The coverage threshold:'))
cdr3 = input('If only pick those antibodies that have peptides in cdr3[Y or N]:')
output_file = input('Need an output file[Y or N]:')

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

with open(pos_file, 'r') as f2:
    pos_dict = {}
    for line in f2:
        line = line.rstrip()
        if line.startswith('>'):
            ID = line.replace('>','')
        else:
            pos_cdr3_st = line.split(' ')[5]
            pos_cdr3_ed = line.split(' ')[6]
            pos_dict[ID] = (pos_cdr3_st,pos_cdr3_ed)

def colored_seq(ref_seq,rare_pep,contigs_hit_pos,cdr3_pos):
    rare_pep_pos = []
    if_rare_cdr3 = False
    for pep in rare_pep:
        start_pos = ref_seq.find(pep)
        end_pos = start_pos + len(pep)
        rare_pep_pos.append((start_pos, end_pos))
    for pos in rare_pep_pos:
        pos_range = range(pos[0],pos[1])
        if cdr3_pos[0] in pos_range or cdr3_pos[1] in pos_range:
            if_rare_cdr3 = True
            hit_range = pos_range
        else:
            hit_range = range(0)
    common_pep_pos = [x for x in contigs_hit_pos if x not in rare_pep_pos]
    colored_seq = ''
    for i in range(len(ref_seq)):
        if any(start <= i < end for start,end in rare_pep_pos):
            colored_seq += red(ref_seq[i])
        elif any(start <= i < end for start,end in common_pep_pos):
            colored_seq += green(ref_seq[i])
        else:
            colored_seq += ref_seq[i]
    return colored_seq,if_rare_cdr3,hit_range

with open(protein_result, 'r') as f3:
    info_dict = {}
    for line in f3:
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

def print_pro(protein):
    pro_style = blue(protein, style='bold+underline')
    ref_seq = info_dict[protein][1][0]
    rare_pep = info_dict[protein][1][1]
    contigs_hit_pos = info_dict[protein][1][2]
    cdr3_pos = (int(pos_dict[protein][0]),int(pos_dict[protein][1]))
    cov = info_dict[protein][0]
    if cov > cov_thres:
        if cdr3 == 'Y':
            if colored_seq(ref_seq,rare_pep,contigs_hit_pos,cdr3_pos)[1]:
                hit_range = colored_seq(ref_seq,rare_pep,contigs_hit_pos,cdr3_pos)[2]
                set_hit = set(hit_range) & set(range(cdr3_pos[0],cdr3_pos[1]))
                if len(set_hit) > 5:
                    print('Protein:'+ pro_style)
                    print('Ref_seq:' + colored_seq(ref_seq,rare_pep,contigs_hit_pos,cdr3_pos)[0])
                    print('Hit_pos:' + str(contigs_hit_pos))
                    print('Coverage:' + str(cov*100) + '%')
                    print(info_dict[protein][2]+'\n')
        else:
            print('Protein:'+ pro_style)
            print('Ref_seq:' + colored_seq(ref_seq,rare_pep,contigs_hit_pos,cdr3_pos)[0])
            print('Hit_pos:' + str(contigs_hit_pos))
            print('Coverage:' + str(cov*100) + '%')
            print(info_dict[protein][2]+'\n')

for protein in info_dict:
    print_pro(protein)
