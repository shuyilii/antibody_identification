#!/usr/bin/python3
import sys
import matplotlib.pyplot as plt

database = sys.argv[1]
analysis_result = sys.argv[2]

with open(database, 'r') as fin:
    size_dict = {}
    for line in fin:
        line = line.rstrip()
        if line.startswith('>'):
            size = line.split(';')[1]
            ID = line.replace('>','').split('|')[0] + '_' + line.replace('>','').split('|')[1]
            size_dict[ID] = size

with open(analysis_result, 'r') as fin:
    read_num_list = []
    for line in fin:
        line = line.rstrip()
        if line.startswith('Protein'):
            ID = line.replace('Protein:', '')
            read_num = size_dict[ID].replace('size=','')
            print(line)
            print('tran_seq_reads_num:' + read_num)
            read_num_list.append(int(read_num))
        else:
            print(line)
count = [[x,read_num_list.count(x)] for x in set(read_num_list)]
pos = [item[0] for item in count]
num = [item[1] for item in count]
count.sort(key = lambda x: x[0])
print('\nreads_num\tantibodies_num\n')
for each in count:
    print(str(each[0]) + '\t' + str(each[1]))
plt.bar(pos, num, width = 1)
plt.xlabel('NGS reads number')
plt.ylabel('antibodies number')
plt.title(analysis_result.replace('.txt',''))
plt.scatter(pos,num,s = 8)
plt.show()
plt.savefig(analysis_result.replace('.txt',''))
