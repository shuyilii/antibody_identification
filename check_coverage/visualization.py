#!/usr/bin/python3


import matplotlib.pyplot as plt

def visualization(each_mapping_list,ref_database_dict):
    max_length = len(max(list(ref_database_dict.values()),key = len))
    reads_num = len(each_mapping_list)
    ax = plt.axes(xlim = (0,max_length),ylim = (0,reads_num + 2))
    i = 1
    for mapping_pos in each_mapping_list:
        ax.arrow(mapping_pos[0], i, mapping_pos[1]-mapping_pos[0], 0, head_width=0.05, head_length=0.1, fc='k', ec='k')
        i += 1
    plt.show()
