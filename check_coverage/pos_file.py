#!/usr/bin/python3

import sys
raw_ngs = sys.argv[1]
full_length_seq = sys.argv[2]

with open(full_length_seq, 'r') as full_length_seq:
    seq_dict = {}
    for line in full_length_seq:
        line = line.rstrip()
        if line.startswith('>'):
            ID = line.replace('>','')
        else:
            seq_dict[ID] = line

with open(raw_ngs, 'r') as raw_ngs:
    next(raw_ngs)
    for line in raw_ngs:
        element_list = line.rstrip().split('\t')
        fr1 = element_list[90].replace(' ','')
        cdr1 = element_list[91].replace(' ','')
        fr2 = element_list[92].replace(' ','')
        cdr2 = element_list[93].replace(' ','')
        fr3 = element_list[94].replace(' ','')
        cdr3 = element_list[95].replace(' ','')
        ID = element_list[27]+ '_' + element_list[34]
        if ID in seq_dict.keys():
            fr1_pos = seq_dict[ID].find(fr1) +1
            cdr1_pos = seq_dict[ID].find(cdr1) +1
            fr2_pos = seq_dict[ID].find(fr2) +1
            cdr2_pos = seq_dict[ID].find(cdr2) +1
            fr3_pos = seq_dict[ID].find(fr3) +1
            cdr3_pos = seq_dict[ID].find(cdr3) +1
            cdr3_end = cdr3_pos + len(cdr3)
            print('>' + ID)
            print(fr1_pos,cdr1_pos,fr2_pos,cdr2_pos,fr3_pos,cdr3_pos,cdr3_end)
