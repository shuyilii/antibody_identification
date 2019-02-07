#!/bin/sh

echo '[This script is for filtering potentially interesting antibodies by getting the antibodies that only in the Php1 or Php2 purification but not in IgG and having at least one distinguishable identified peptide]'
echo 'input file1:'
read file1
echo 'input file2:'
read file2
echo -n 'input filter number for distinguishable identified peptide:'
read fil_num

./uniq_pep.py $file1 > ${file1%filter.txt}uniq_pep.txt

./filter.py $file1 ${file1%filter.txt}uniq_pep.txt $fil_num > ${file1%filter.txt}mark.txt

./filter_1.py $file1 $file2 > ${file1%filter.txt}only.txt

./filter_2.py ${file1%filter.txt}mark.txt ${file1%filter.txt}only.txt > ${file1%filter.txt}filter_result.txt

echo 'result file:' ${file1%filter.txt}filter_result.txt 'is in your data dict'
