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
    parser.add_argument('--siRNA-file', help='File path to siRNA data (CSV)',
                        type=str)
    parser.add_argument('--gene-file', help='File path to gene list (CSV)',
                        type=str)
    parser.add_argument('--on-target-effect',
                        help='Target prediction for on-target effects',
                        type=float, default=0.75)
    parser.add_argument('--outdir', help='Output directory filepath',
                        type=str)

    args = parser.parse_args()
    missing = (args.cps_dir is None or
               args.siRNA_file is None or
               args.gene_file is None or
               args.outdir is None)
    if missing:
        raise ValueError("Missing arguments")

    genes = pd.DataFrame.from_csv(args.gene_file, index_col=False)
    genes['GeneSynonymsSplit'] = (genes
                                  .GeneSynonyms
                                  .apply(lambda x: x.split(';')))
    genes.index.name = 'SparseGeneID'
    gene_dict = dict(zip(genes.GeneID, genes.index))

    siRNAs = pd.DataFrame.from_csv(args.siRNA_file, index_col=False)
    siRNAs.index.name = 'SparseSiRNAID'
    ln = len(gene_dict)
    for geneID, geneName in zip(siRNAs.geneID, siRNAs.gene):
        if geneID not in gene_dict:
            gene_dict[geneID] = ln
            genes = pd.concat([genes, pd.DataFrame({'GeneID': int(geneID),
                                                    'GeneName': geneName,
                                                    'GeneSynonyms': '-',
                                                    'GeneSynonymsSplit':
                                                    [['-']]},
                                                   index=[int(ln)])])
            ln += 1
    gene_dict = dict(zip(genes.GeneID, genes.index))

    sparse_mat = dok_matrix((len(siRNAs), len(genes)), dtype=np.float32)
    for si, si_id, gn_id in zip(siRNAs.siRNAID, siRNAs.index, siRNAs.geneID):
        siRNA_fname = os.path.join(args.cps_dir,
                                   '%s.tsv' % si)
        di = sparsify_cps(siRNA_fname,
                          si_id,
                          gene_dict)
        for k, v in di.iteritems():
            sparse_mat[k[0], k[1]] = max(1. - (2 ** v), 0.0)
        sparse_mat[si_id, gene_dict[gn_id]] = 0.75
    csr_mat = sparse_mat.tocsr()

    sparse_outfile = os.path.join(args.outdir, 'X.csr')
    genes_outfile = os.path.join(args.outdir, 'genes.csv')
    siRNAs_outfile = os.path.join(args.outdir, 'siRNAs.csv')

    save_sparse_csr(sparse_outfile, csr_mat)
    genes.to_csv(genes_outfile, index=True)
    siRNAs.to_csv(siRNAs_outfile, index=True)
