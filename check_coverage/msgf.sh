#!/bin/sh


for file in ../monoclonals/*.mgf
do
  java -jar ../software/msgf+/MSGFPlus.jar -s $file -d ../monoclonals/monoclonal.fasta -t 10ppm -tda 1 -o ${file%.mgf}.mzid -m 3 -n 10 -thread 8
  java -Xmx3500M -cp ../software/msgf+/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i ${file%.mgf}.mzid -unroll 1
done
