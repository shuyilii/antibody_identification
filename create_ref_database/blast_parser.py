#!/usr/bin/python3

'''
The function is used for parsering the blast result of the trimmed sequence to the IMBT constant region database.
The input is the blast outcome, the output is a dictionary {query_id : [constant_id, position]}
'''

def blast_parser(blast_outcome):

    query_list = []
    results_list = []
    pos_constant = []

    with open(blast_outcome,'r') as fin:
        for line in fin:
            line = line.rstrip()

            if line.startswith('Query'):
                element_q = line.split('  ')
                if len(element_q) < 2:
                    query = element_q[0].replace('Query= ','').split(';')[0]
                    query_list.append(query)
                    #print(query)

            if line.startswith('*'):
                results_list.append('null')
                pos_constant.append('null')
                #print('0'+'\n'+'0')

            if line.startswith('>'):
                constant_id = line.replace('> ','')
                results_list.append(constant_id)
                Check_multi_hits_sbjct = False
                #print(constant_id)

            if line.startswith(' Identities'):
                power = line.split(' ')
                match_base_num = power[3].split('/')[1]

            if line.startswith('Sbjct'):
                element_c = line.split('  ')
                if match_base_num == '60':
                    if len(element_c) == 4:
                        pos_constant.append(element_c[3])
                        #print(element_c[3])
                    else:
                        pos_constant.append(element_c[4])
                        #print(element_c[4])
                else:
                    if len(element_c) == 4:
                        pos_constant.append(element_c[3])
                        Check_multi_hits_sbjct = True
                        #print(element_c[3])
                    else:
                        if int(element_c[4]) < int(element_c[1]) + 59:
                            if Check_multi_hits_sbjct == False:
                                pos_constant.append(element_c[4])
                                #print(element_c[4])
    # print(len(query_list))
    # print(len(results_list))
    # print(len(pos_constant))

    temp_list = []
    for constant_id, pos in zip(results_list, pos_constant):
        temp_list.append([constant_id, pos])
    results_dict = dict(zip(query_list, temp_list))

    return results_dict

#
# import sys
# blast_outcome = sys.argv[1]
# blast_parser(blast_outcome)
