# Predict siRNA-to-gene off-target relation matrices using TargetScan

Authors: Fabian Schmich (fabian.schmich@bsse.ethz.ch) and Mason Victors (mason.victors@recursionpharma.com)

The following steps will guide you through a set scripts used for the prediction of siRNA off-targets and the construction of siRNA-to-gene target relation matrices as used for input to the gespeR model (www.cbg.ethz.ch/software/gespeR). It is required to install TargetScan and TargetScan context+ (http://www.targetscan.org/). Please adjust paths in each script according to your installation.

1. Prepare all seed sequences running the script "make_seeds.py".

   This script creates one file for each n-mer seed sequence, suitable for input to TargetScan's targetscan_60.pl input script

2. Predict all seed matches running script "targetscan_60.py".

   This script loops over all seeds created in step 1. For consecutive predictions, i.e. for multiple libraries, step 1 and 2 only need to be done once. We have provided a file containing the 3' UTR sequences with with removed gaps for hg19 in the ./data folder.

3. Predict the context+ scores for all matches found in step 2 running script "run_predicttargets.py".

   This script should take your results from step 2 and calculate the context+ score for each gene for each file produced by step 2 (corresponding to the seeds).

4. Construct the siRNA-to-gene target relation matrix from all context+ score predictions in step 3 running script "make_X.sh".

   This script creates a percent knockdown csv for all seeds generated and all genes, appending any missing critical genes.


Example:

```mkdir data/seeds
mkdir data/tscanout
mkdir data/cpscores
python make_seeds.py data/seeds
python targetscan_60.py --utr-seq-file data/hg19_ucsc_3p.tsv --seeds-directory data/seeds --out-directory data/tscanout --n-workers 10
python run_predicttargets.py --tscan-outdir data/tscanout --utr-file data/hg19_ucsc_3p.tsv --ta-sps-file data/TA_SPS_by_seed_region.txt --ref-seq-file data/refseq.tsv --n-workers 10 --outdir data/cpscores
python make_X.py --cps-dir data/cpscores --outdir data --critical-genes data/critical_genes.csv
