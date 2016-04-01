#!/bin/bash

#sudo cpan App::cpanminus
#. ~/.bashrc
#sudo cpanm Statistics::Lite
#sudo cpanm Bio::TreeIO
#
#sudo apt-add-repository ppa:j-4/vienna-rna
#sudo apt-get update
#sudo apt-get install vienna-rna

#mkdir -p "$1/TargetScan7"
#wget --directory-prefix="$1/TargetScan7" http://www.targetscan.org/vert_70/vert_70_data_download/targetscan_70.zip
#cd "$1/TargetScan7"
#unzip targetscan_70.zip
#rm targetscan_70.zip
echo "export PATH=$1/TargetScan7:\$PATH" >> ~/.bashrc
#cd ~
#
#wget --directory-prefix="$1" http://www.targetscan.org/vert_70/vert_70_data_download/targetscan_70_BL_PCT.zip
#cd "$1"
#unzip targetscan_70_BL_PCT.zip
#rm targetscan_70_BL_PCT.zip
echo "export PATH=$1/TargetScan7_BL_PCT:\$PATH" >> ~/.bashrc
#cd ~
#
#wget --directory-prefix="$1" http://www.targetscan.org/vert_70/vert_70_data_download/TargetScan7_context_scores.zip
#cd "$1"
#unzip TargetScan7_context_scores.zip
#rm TargetScan7_context_scores.zip
echo "export PATH=$1/TargetScan7_context_scores:\$PATH" >> ~/.bashrc
#cd ~
