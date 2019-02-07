#!/usr/bin/python3

import sys
from itertools import islice

Php2_mark = sys.argv[1]
Php2_only = sys.argv[2]
with open(Php2_only, 'r') as f1:
    protein_list = []
    for line in f1:
        line = line.rstrip()
        if line.startswith('Protein'):
            pro = line.replace('Protein:','')
            protein_list.append(pro)
with open(Php2_mark, 'r') as f2:
    for line in f2:
        line = line.rstrip()
        if line.startswith('Protein'):
            pro = line.replace('Protein:','')
            if pro in protein_list:
                print('Protein:' + pro)
                print(''.join(islice(f2, 7)))
