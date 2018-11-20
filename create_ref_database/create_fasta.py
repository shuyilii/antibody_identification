#!/usr/bin/python3
'''
This is for create a fasta file from the csv file containing information:
[isosubtype, part_igh_id, trimmed_sequence, v_segment,stop_codon,productive, v_start pre_seq_aa_q,fr1_seq_aa_q,cdr1_seq_aa_q,fr2_seq_aa_q,cdr2_seq_aa_q,fr3_seq_aa_q,cdr3_seq_aa_q,post_seq_aa_q]
Some IMGT v_segment access number have been updated in the database, so in this script it is updated accordingly.
input: csv file / output: fasta file
> ID | v_segment | isotype | v_start | q_start
nt_seq OR aa_seq
'''

import sys

mode = sys.argv[3]
database_trans_dic = {'c':'38-4','f':'69-2','d':'38-3','h':'69-1','b':'38-2','a':'10-1'}
with open(sys.argv[1],'r') as fin, open(sys.argv[2],'w') as fout:
    for line in fin:
        line = line.rstrip()
        elements_list = line.split('\t')
        isotype = elements_list[0]
        ID = elements_list[1]
        nt_seq = elements_list[2]
        v_segment = elements_list[3]
        if v_segment[6] in list(database_trans_dic.keys()):
            v_segment = v_segment[0:6] + database_trans_dic[v_segment[6]] + v_segment[7:]
        if v_segment == 'IGHV2-5*10':
            v_segment = 'IGHV2-5*02'
        if v_segment == 'IGHV2-5*07':
            v_segment = 'IGHV2-5*04'
        if v_segment == 'IGHV3-43D*01':
            v_segment = 'IGHV3-43D*03'
        stop_codon = elements_list[4]
        productive = elements_list[5]
        q_start = elements_list[6]
        v_start = elements_list[7]
        aa_seq = ''.join(elements_list[8:]).replace(' ', '')
###unprodictive sequences and the sequences that contains stop_codon are filtered. And we only takes IGHG isotype into account
        if productive == 't' and stop_codon == 'f':
            if aa_seq.find('*') < 0:
                if isotype.find('IGHG') >= 0:
                    if mode == 'nt':
                        print('>' + ID + '|' + v_segment + '|' + isotype + '|' + v_start +'|' + q_start + '\n' + nt_seq, file = fout )
                    if mode == 'aa':
                        print('>' + ID + '|' + v_segment + '|' + isotype + '\n' + aa_seq, file = fout)
