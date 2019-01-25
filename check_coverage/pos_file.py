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

database_trans_dic = {'c':'38-4','f':'69-2','d':'38-3','h':'69-1','b':'38-2','a':'10-1'}
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
        v_segment = element_list[34]
        if v_segment[6] in list(database_trans_dic.keys()):
            v_segment = v_segment[0:6] + database_trans_dic[v_segment[6]] + v_segment[7:]
        if v_segment == 'IGHV2-5*10':
            v_segment = 'IGHV2-5*02'
        if v_segment == 'IGHV2-5*07':
            v_segment = 'IGHV2-5*04'
        if v_segment == 'IGHV3-43D*01':
            v_segment = 'IGHV3-43D*03'
        ID = element_list[27]+ '_' + v_segment
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
