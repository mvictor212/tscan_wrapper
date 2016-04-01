#!/usr/bin/env python

# Create all seeds of length SEED_LEN. These files will serve as input to TargetScan's base script (targetscan_70.pl).
#
# Authors: Fabian Schmich (fabian.schmich@bsse.ethz.ch) and Mason Victors
# (mason.victors@recursionpharma.com)
#

# Libraries
import os
from itertools import product


def make_seeds(out_directory, species, seed_length):
    alphabet = ["A", "C", "G", "U"]
    seeds = [''.join(i) for i in product(alphabet, repeat=seed_length)]
    for s in seeds:
        f = open("%s/%s.txt" % (out_directory, s), "w")
        f.write("%s\t%s\t%s\n" % (s, s, species))
        f.close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--species', help='Species of interest',
                        type=str, default='9606')
    parser.add_argument('--seed-length', help='Length of seeds',
                        type=int, default=7)
    parser.add_argument('directory', help='Directory to save seed sequence files')

    args = parser.parse_args()
    make_seeds(args.directory, args.species, args.seed_length)
