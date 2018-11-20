#!/bin/sh

echo 'create a full length aa/nt antibody sequences(with both vdj region and CH1 constant region) from the trimmed sequence'
echo 'input file:'
read file

##extract sequences
cat $file | cut -f 27,28,32,35,42,44,50,52,90-97 > ${file}_temp

##extract the trimmed sequence
../antibody_identification/create_ref_database/create_fasta.py ${file}_temp nt_seq.fasta nt

##start blast for the c_region (Nucleotide-Nucleotide Version: BLAST 2.7.1+)
printf '\nStart blast'
mkdir ../blast_database
cd ../blast_database
makeblastdb -in ../database/Ig_c.fasta -out Ig_c_Nt -dbtype nucl
mkdir ../blast
cd ../blast
blastn -query ../NGS/nt_seq.fasta -db ../blast_database/Ig_c_Nt -evalue 1e-10 -out blast_results.blastn -num_alignments 1 -num_descriptions 1 -num_threads 4 -gapopen 99 -gapextend 9
printf '\nEnd blast\n'

## add the missing c_region and v_region to the trimmed_sequence
printf '\nStart creating database\n'
../antibody_identification/create_ref_database/fulllength.py ../database/Ig_v.fasta ../database/Ig_c.fasta blast_results.blastn ../NGS/nt_seq.fasta > full_length.fasta

##translated nucleotide sequence to the amino acid sequence (Version: EMBOSS:6.6.0.0)
transeq -frame 1 -sformat pearson -sequence full_length.fasta -outseq trans_fulllength_seq.fasta 2> err.log

##sequence dereplication (Version:usearch v11.0.667_i86linux32)
##Sequences are compared letter-by-letter and must be identical over the full length of both sequences (substrings do not match).
usearch -fastx_uniques trans_fulllength_seq.fasta -sizeout -fastaout trans_seq.derep.fasta 2>> err.log

##add other database to control FDR (ref:http://www.peptideatlas.org/thisp/    level2 database)
cat trans_seq.derep.fasta ../database/PA_THISP_Level2_2018-11-01.fasta > ref_database_temp.fasta

##transform the database to the one that can be processed by the search engine
../antibody_identification/create_ref_database/trans_format_database.py ref_database_temp.fasta > ref_database.fasta
printf '\nEnd creating database\n'
