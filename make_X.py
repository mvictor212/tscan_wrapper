#!/usr/bin/python

# Run makeX.py on the output directory of predicted context+ scores
# for siRNA off-targets.
#
# Authors: Fabian Schmich (fabian.schmich@bsse.ethz.ch) and Mason Victors
# (mason.victors@recursionpharma.com)
#

import os
import numpy as np
import pandas as pd
from scipy.sparse import dok_matrix, csr_matrix


def save_sparse_csr(fname, arr):
    np.savez(fname, data=arr.data, indices=arr.indices,
             indptr=arr.indptr, shape=arr.shape)


def load_sparse_csr(fname):
    loader = np.load(fname)
    return csr_matrix((loader['data'], loader['indices'], loader['indptr']),
                      shape=loader['shape'])


def sparsify_cps(cpsfile, si_id, gene_dict):
    df = pd.DataFrame.from_csv(cpsfile, sep='\t', index_col=False)
    return {(si_id, gene_dict[gn]): score
            for gn, score in zip(df.GeneID, df.CPSmean)}


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--cps-dir',
                        help='Context+ Score Predictions directory',
                        type=str)
    parser.add_argument('--outdir', help='Output directory filepath',
                        type=str)
    parser.add_argument('--critical-genes', help='Path to csv of critical entrez gene ids',
                        type=str)

    args = parser.parse_args()
    missing = (args.cps_dir is None or
               args.outdir is None or
               args.critical_genes is None)
    if missing:
        raise ValueError("Missing arguments")

    critical_genes = pd.DataFrame.from_csv(args.critical_genes, index_col=None)
    target_relations_dict = {}
    for i, seq in enumerate(os.listdir(args.cps_dir)):
        if (i % 1000) == 0:
            print("Iteration: %d" % i)
        df = pd.DataFrame.from_csv(os.path.join(args.cps_dir, seq),
                                   sep='\t', index_col=0)
        try:
            target_relations_dict[seq.split('.')[0]] = df.CPSmean.apply(lambda x: max(1 - 2 ** x, 0.))
        except:
            print(df.columns)
            print(seq)

    target_relations_df = pd.DataFrame(target_relations_dict).T
    target_relations_df.fillna(0., inplace=True)
    target_relations_df.index.name = 'entrez_gene_id'
    target_relations_df.to_csv(os.path.join(args.outdir, 'seed_by_gene.csv'))
    missing_genes = list(set(critical_genes
                             .entrez_gene_id
                             .values
                             .astype(int))
                         .difference(set(target_relations_df
                                         .columns
                                         .astype(int))))
    missing_genes_df = pd.DataFrame(0, columns=missing_genes,
                                    index=target_relations_df.index)
    target_relations_df = pd.concat([target_relations_df, missing_genes_df], axis=1)
    target_relations_df.to_csv(os.path.join(args.outdir, 'seed_by_gene_extended.csv'))
