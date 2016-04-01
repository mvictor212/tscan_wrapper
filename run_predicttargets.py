#!/usr/bin/env python

# Authors: Fabian Schmich (fabian.schmich@bsse.ethz.ch) and Mason Victors
# (mason.victors@recursionpharma.com)
#
import os
import pandas as pd
import multiprocessing as mp

import predicttargets as pt


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--siRNA-file', help='File path to siRNA data',
                        type=str)
    parser.add_argument('--tscan-outdir',
                        help='targetscan_60.pl output filepath',
                        type=str)
    parser.add_argument('--utr-file', help='UTR Seq filepath',
                        type=str)
    parser.add_argument('--ref-seq-file',
                        help='Reference file for transcript to gene lookup',
                        type=str)
    parser.add_argument('--n-workers', help='Number of parallel processes',
                        type=int, default=12)
    parser.add_argument('--outdir', help='Output directory filepath',
                        type=str)

    args = parser.parse_args()
    missing = (args.siRNA_file is None or
               args.tscan_outdir is None or
               args.utr_file is None or
               args.ref_seq_file is None or
               args.outdir is None)
    if missing:
        raise ValueError("Missing arguments")

    translation_dict = pt.get_translation_dict(args.ref_seq_file)
    df = pd.DataFrame.from_csv(args.siRNA_file, index_col=False)

    def foo(tup):
        pt.write_target_frame(tup[0], tup[1],
                              args.tscan_outdir,
                              args.utr_file,
                              args.outdir,
                              translation_dict=translation_dict)

    genes = pd.DataFrame(columns=['GeneID', 'GeneName', 'GeneSynonyms'])
    if args.n_workers == 1:
        for tup in zip(df['siRNAID'], df['sequence']):
            genes = pd.concat([genes, foo(tup)]).drop_duplicates()
    else:
        p = mp.Pool(args.n_workers)
        ret = p.map(foo, zip(df['siRNAID'], df['sequence']))
        p.terminate()
        genes = pd.concat(ret).drop_duplicates()
    genes.to_csv(os.path.join(args.outdir,
                              'genes.csv'),
                 index=False)
