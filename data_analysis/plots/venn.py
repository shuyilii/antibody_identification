#!/usr/bin/python3

import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles

A_filter = sys.argv[1]
B_filter = sys.argv[2]
C_filter = sys.argv[3]

def get_pro(file):
    with open(file, 'r') as file:
        pro_fil = []
        for line in file:
            line = line.rstrip()
            if line.startswith('Protein'):
                pro = line.replace('Protein:','')
                pro_fil.append(pro)
    return pro_fil
A = set(get_pro(A_filter))
B = set(get_pro(B_filter))
C = set(get_pro(C_filter))
AB = A & B
AC = A & C
BC = B & C
AorB = A|B
AorC = A|C
BorC = B|C
ABC = AB & C
all = AorB | C
subset = (len(all-BorC),len(all-AorC),len(AB-C),len(all-AorB),len(AC-B),len(BC-A),len(ABC))

plt.figure(figsize = (10,10))
v = venn3(subsets = subset, set_labels = ('Php1', 'Php2', 'IgG'))
c = venn3_circles(subsets = subset, linestyle='dashed')
plt.title("Venn diagram of identified antibodies")
plt.show()
