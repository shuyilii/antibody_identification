#!/usr/bin/python3
import sys


with open(sys.argv[1],'r') as fin:
    for line in fin:
        line = line.rstrip()
        if line.startswith('>'):
            elements_list = line.split('|')
            ID = '_'.join(elements_list[0:2])
            print(ID)
        else:
            line = line.replace('*','X')
            print(line)
