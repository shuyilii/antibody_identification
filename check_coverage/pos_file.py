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

database_trans_dic = {'c':'38-4','f':'69-2','d':'38-3','h':'69-1','b':'38-2','a':'10-1',
                      'IGHV2-5*10':'IGHV2-5*02','IGHV2-5*07':'IGHV2-5*04','IGHV3-43D*01':'IGHV3-43D*03'}
with open(raw_ngs, 'r') as raw_ngs:
    next(raw_ngs)
    for line in raw_ngs:
        element_list = line.rstrip().split('\t')
        v_segment = element_list[34]
        if v_segment[6] in database_trans_dic.keys():
            v_segment = v_segment[0:6] + database_trans_dic[v_segment[6]] + v_segment[7:]
        if v_segment in database_trans_dic.keys():
            v_segment = database_trans_dic[v_segment]
        ID = element_list[27]+ '_' + v_segment
        region_list = []
        for a in range(90,96):
            region_seq = element_list[a].replace(' ','')
            region_list.append(region_seq)
        if ID in seq_dict.keys():
            pos_list = []
            for region_seq in region_list:
                pos_list.append(seq_dict[ID].find(region_seq) +1)
            cdr3_end = pos_list[-1] + len(region_list[-1])
            pos_list.append(cdr3_end)
            print('>' + ID)
            print(pos_list)
