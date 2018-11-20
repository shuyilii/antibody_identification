#!/usr/bin/python3
#
import sys

####transformat database ig_c :transform the ID line
with open(sys.argv[1],'r') as fin:
    for line in fin:
        line = line.rstrip()
        if line.startswith('>'):
            elements_list = line.split('|')
            ID = '_'.join(elements_list[0:5])
            print(ID)
        else:
            print(line)


#####check blast_parser
# def RepresentsInt(s):
#     try:
#         int(s)
#         return True
#     except ValueError:
#         return False
# #
# i = 0
# with open(sys.argv[1],'r') as fin:
#     for line in fin:
#         line = line.rstrip()
#         i += 1
#         if i%4 == 1:
#             tag = line
#         if i%4 == 0:
#             if RepresentsInt(line):
#                 pass
#             else:
#                 print(tag)
#                 sys.exit()
#
# with open(sys.argv[1],'r') as fin:
#     for line in fin:
#         line = line.rstrip()
#         if RepresentsInt(line):
#             if int(line) <= 60 and int(line) != 0:
#                 print(line)


###check if the translated seq is identical as the original one
# dict_1 = {}
# dict_2 = {}
# with open(sys.argv[1], 'r') as f1:
#     for line in f1:
#         line = line.rstrip()
#         if line.startswith('>'):
#             ID = line.replace('>','').split('|')[0]
#         else:
#             dict_1[ID] = line
# with open(sys.argv[2], 'r') as f2:
#     for line in f2:
#         line = line.rstrip()
#         if line.startswith('>'):
#             ID = line.replace('>','').split('|')[0]
#         else:
#             dict_2[ID] = line
#
# for each in dict_2:
#     if dict_2[each].find(dict_1[each]) < 0:
#         print(each)

###transformat translated_seq : transform '*'-'Z'
# with open(sys.argv[1], 'r') as f1:
#     for line in f1:
#         line = line.rstrip()
#         if line.startswith('>'):
#             print(line)
#         else:
#             line = line.replace('*','Z')
#             print(line)
