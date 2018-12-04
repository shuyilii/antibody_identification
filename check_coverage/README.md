* mapping the identified peptides to the reference protein sequence

* check the mapping coverage

* visualize the mapping process

usage:

run_check_coverage.py [-h help] [-r ref] [-i input] [-c1 cutoff1] [-c2 cutoff2]

sample usage:

run_check_coverage.py -r ref_database.fasta -i data1.tsv data2.tsv data3.tsv (multiple files here) -c1 0.01 -c2 0.18

use 'run_check_coverage.py -h' for more information
