#!/usr/bin/python3

import sys
#import pandas as pd

tsv = sys.argv[1]
mgf = sys.argv[2]

# tsv_file = pd.read_table(tsv)
# title_list = list(tsv_file['Title'])

with open(mgf, 'r') as mgf:
    intensity_dict = {}
    for line in mgf:
        line = line.rstrip()
        if line.startswith('TITLE='):
            title = line.replace('TITLE=','')
        if line.startswith('PEPMASS='):
            prec_intensity = line.split(' ')[1]
            intensity_dict[title] = prec_intensity
with open(tsv, 'r') as tsv:
    print('SpecFile\tSpecID\tScanNum\tTitle\tFragMethod\tPrecursor\tIsotopeError\tPrecursorError(ppm)\tCharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\tSpecEValue\tEValue\tQValue\tPepQValue\tPrecursor intensity')
    next(tsv)
    for line in tsv:
        line = line.rstrip()
        title = line.split('\t')[3]
        print(line + '\t' + intensity_dict[title])
