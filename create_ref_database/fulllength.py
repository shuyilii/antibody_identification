#!/usr/bin/python3

import sys
from blast_parser import blast_parser

ig_v = sys.argv[1]
ig_c = sys.argv[2]
blast_results = sys.argv[3]
nt_fasta = sys.argv[4]

with open(ig_v,'r') as ig_v, open(nt_fasta, 'r') as sample:
    ref = {}
    full_v_dict = {}
    for line in ig_v:
        line = line.rstrip()
        if line.startswith('>'):
            ref_elements_list = line.replace('>','').split('|')
            Igv_ref = ref_elements_list[1]
        else:
            ref_seq = line
            ref[Igv_ref] = ref_seq

    for line in sample:
        line = line.rstrip()
        if line.startswith('>'):
            smp_elements_list = line.replace('>','').split('|')
            Igv_smp = smp_elements_list[1]
            v_start = int(smp_elements_list[3])
            q_start = int(smp_elements_list[4])
            if Igv_smp in list(ref.keys()):
                if v_start-q_start > 0:
                    if q_start < 20:
                        v_remain_seq = ref[Igv_smp][0:v_start-q_start]
                        smp_id = line.replace('>','')
                        chimera_seq = False
                    else:
                        v_remain_seq = ref[Igv_smp][0:v_start-1]
                        chimera_seq = True
                else:
                    v_remain_seq = ''
                    smp_id = line.replace('>','')
            else:
                print ('error: Igv not in the database')
                sys.exit()
        else:
            if chimera_seq == False:
                smp_seq = v_remain_seq + line
                full_v_dict[smp_id] = smp_seq.upper()
            else:
                smp_seq = v_remain_seq + line[q_start-1:]
                full_v_dict[smp_id] = smp_seq.upper()


Igc_pos_dict = blast_parser(blast_results)

with open(ig_c,'r') as ig_c:
    ref = {}
    for line in ig_c:
        line = line.rstrip()
        if line.startswith('>'):
            Igc_ref = line.replace('>','')
        else:
            ref_seq = line
            ref[Igc_ref] = ref_seq
for smp_id in full_v_dict.keys():
    if smp_id in Igc_pos_dict.keys():
        try:
            int(Igc_pos_dict[smp_id][1])
            pos_constant = int(Igc_pos_dict[smp_id][1])
            Igc_ref = Igc_pos_dict[smp_id][0]
            c_remain_seq = ref[Igc_ref][pos_constant:]
            fulllength_seq = full_v_dict[smp_id] + c_remain_seq
            print('>' + smp_id)
            print(fulllength_seq.upper())
        except ValueError:
            pass
