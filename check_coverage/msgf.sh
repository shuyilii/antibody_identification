#!/bin/sh


for file in *.mgf
do
  java -jar ../software/msgf+/MSGFPlus.jar -s $file -d ../../database/ref_database.fasta -t 10ppm -tda 1 -inst 3 -e 5 -o ${file%.mgf}.mzid -m 3 -n 10 -thread 8
  java -Xmx3500M -cp ../software/msgf+/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i ${file%.mgf}.mzid -unroll 1
done
