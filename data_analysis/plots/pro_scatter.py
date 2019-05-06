#!/usr/bin/python3

import sys
from matplotlib import pyplot as plt

rank_file = sys.argv[1]

d1 = [[],[]]
d2 = [[],[]]
d3 = [[],[]]
with open(rank_file) as fin:
    next(fin)
    i = 0
    for line in fin:
        element = line.rstrip().split('\t')
        antibodies = element[0]
        Php2 = int(element[1])
        IgG = int(element[2])
        diff = int(element[3])
        i += 1
        if i < 50:
            d1[0].append(Php2)
            d1[1].append(IgG)
        else:
            if diff > 0:
                d2[0].append(Php2)
                d2[1].append(IgG)
            else:
                d3[0].append(Php2)
                d3[1].append(IgG)
data = (d3,d2,d1)
colors = ('k','g','r')
groups = ('not_interested','interested','top50')
for data, color, group in zip(data, colors, groups):
    x, y = data
    plt.scatter(x, y, c = color,label = group,s = 10)
plt.legend(loc = 'upper right')
plt.xlabel('Php2', fontsize = 10)
plt.ylabel('IgG', fontsize = 10)
plt.ylim(0, 4000)
plt.xlim(0, 4000)
plt.title('scatter plot of antibodies intensity')
plt.show()
