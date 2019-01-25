#!/usr/bin/python3
import sys

with open(sys.argv[1],'r') as fin:
    next(fin)
    for line in fin:
        line = line.rstrip()
        elements_list = line.split('\t')
        peptide = ''.join([i for i in elements_list[9] if not i.isdigit()]).replace('+.','')
        protein = elements_list[10]
        Q_value = float(elements_list[15])
        if protein.find('IGHV') >= 0 and Q_value <= 0.5:
            print(peptide + '\t' + protein + '\t' + str(Q_value))
