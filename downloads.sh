#!/bin/bash

wget --output-document="$1/ORF_Sequences.txt.zip" http://www.targetscan.org/vert_70/vert_70_data_download/ORF_Sequences.txt.zip
unzip "$1/ORF_Sequences.txt.zip"

wget --output-document="$1/UTR_Sequences.txt.zip" http://www.targetscan.org/vert_70/vert_70_data_download/UTR_Sequences.txt.zip
unzip "$1/UTR_Sequences.txt.zip"

sed '1,1d' UTR_Sequences.txt | cut -f1,4,5 > UTR_sequences_all.txt
