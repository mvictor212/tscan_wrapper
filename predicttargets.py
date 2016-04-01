#!/usr/bin/env python

# Run TargetScan's context score script (targetscan_60_context_scores.pl) on
# the output of TargetScan's base script (targetscan_60.pl). The script takes
# as input the id and sequence of an siRNA, as well as the directory of seeds
# scanned with TargetScan's base script and a translation file from transcripts
# to genes, and outputs TargetScan context score predictions to an output
# folder. Remember, we cannot use aligned UTRs, since siRNAs are not endogenous
# entities.
#
# Authors: Fabian Schmich (fabian.schmich@bsse.ethz.ch) and Mason Victors
# (mason.victors@recursionpharma.com)
#

# Libraries
import sys
import os.path
import optparse
import tempfile
import subprocess
import pandas as pd
from numpy import mean

# Paths
PROC = subprocess.Popen(['bash',
                         '-c',
                         'which targetscan_60_context_scores.pl'],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                        stdin=subprocess.PIPE)
STDOUT, STDERR = PROC.communicate()
CS_SCRIPT = STDOUT.strip()
CS_SCRIPT_DIR = os.path.dirname(CS_SCRIPT)


class Target:
    "Class to represent a siRNA target gene"

    def __init__(self, id, name=None, synonyms=[],
                 scores=[], transcripts=[]):
        self.id = id
        self.name = name
        self.synonyms = synonyms
        self.scores = scores
        self.transcripts = transcripts

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return self.id == other.id

    def __repr__(self):
        if self.name is not None:
            return "\t".join([str(self.id),
                              self.name,
                              ';'.join(self.synonyms)])
        else:
            return str(self.id)

    def strextd(self):
        return "\t".join([str(self.id),
                          self.name,
                          ";".join(self.synonyms),
                          ";".join(self.transcripts),
                          ";".join([str(e) for e in self.scores]),
                          str(mean(self.scores))])


# Takes an siRNAs and produces the input files required by context_score
# script, "miRNA file" and "PredictedTargets" file by matching the sRNA
# seed to pre-scanned seeds.
def prepare(sirna_id, sirna_seq, targetscan_60_outdir, species="9606"):
    # Mature sequence file
    mature_sequence_file = tempfile.NamedTemporaryFile("w",
                                                       prefix="seq_",
                                                       delete=False)
    mature_sequence_file.write("\t".join(["miRNA_family_ID",
                                          "Species_ID",
                                          "MiRBase_ID",
                                          "Mature_sequence"]) + "\n")
    mature_sequence_file.write("\t".join([sirna_id,
                                          str(species),
                                          sirna_id,
                                          sirna_seq]) + "\n")
    mature_sequence_file.close()
    # Seed predictions file
    ts_predictions_file = tempfile.NamedTemporaryFile("w",
                                                      prefix="targets_",
                                                      delete=False)
    seed_filename = os.path.join(targetscan_60_outdir,
                                 'tscan.%s.tsv' % sirna_seq[1:8])
    print sirna_seq[1:8]
    with open(seed_filename, "r") as seedsf:
        ts_predictions_file.write(seedsf.readline())  # header
        for line in seedsf:
            lsp = line.split("\t")
            lsp[1] = sirna_id
            ts_predictions_file.write("\t".join(lsp))
    ts_predictions_file.close()
    return(mature_sequence_file.name, ts_predictions_file.name)


def predict(mat_seq_file_name, ts_pred_file_name, utr_file):
    contextplus_score_file = tempfile.NamedTemporaryFile("w",
                                                         prefix="cps_",
                                                         delete=True)
    cps_fname = contextplus_score_file.name
    contextplus_score_file.close()
    olddir = os.getcwd()
    os.chdir(CS_SCRIPT_DIR)
    print mat_seq_file_name, ts_pred_file_name, cps_fname, utr_file, CS_SCRIPT
    from subprocess import call
    call([CS_SCRIPT, mat_seq_file_name, utr_file,
          ts_pred_file_name, cps_fname])
    os.chdir(olddir)
    return (cps_fname)


def get_translation_dict(ref_seq_file):
    with open(ref_seq_file, "r") as translf:
        tr = dict()
        for line in translf:
            lsp = line.split("\t")
            lsp = [el.strip() for el in lsp]
            if lsp[0] not in tr:
                tr[lsp[0]] = [[lsp[1], lsp[2], lsp[3]]]
            else:
                tr[lsp[0]].append([lsp[1], lsp[2], lsp[3]])
    return tr


def get_transcript_dict(csfile):
    df = pd.DataFrame.from_csv(csfile, sep='\t', index_col=False)
    transcripts = {k: v['context+ score'].values.tolist()
                   for k, v in df.groupby('Gene ID')}
    return transcripts


def get_targets(transcript_dict, translation_dict):
    # Translate
    targets = []
    for trkey in transcript_dict.keys():
        trname = trkey
        if "." in trname:
            trname = trname[:trname.index(".")]  # Strip everything after "."
        if trname in translation_dict:
            for [gid, gn, syn] in translation_dict[trname]:
                curgene = Target(id=gid, name=gn, synonyms=syn.split('|'))
                if curgene in targets:
                    g_ind = targets.index(curgene)
                    targets[g_ind].scores.append(sum(transcript_dict[trkey]))
                    targets[g_ind].transcripts.append(trname)
                else:
                    curgene.scores = [sum(transcript_dict[trkey])]
                    curgene.transcripts = [trname]
                    targets.append(curgene)
        else:
            print("WARN: cannot translate transcript %s" % trname)
    target_frame = pd.DataFrame([target.strextd().split('\t')
                                 for target in targets],
                                columns=['GeneID', 'GeneName',
                                         'GeneSynonyms',
                                         'Transcripts', 'CPS',
                                         'CPSmean'])
    return target_frame


def write_target_frame(si_id, seq, tscan_outdir, utr_file,
                       outdir, ref_seq_file=None,
                       translation_dict=None):
    print si_id, seq
    mat_seq_file, ts_pred_file = prepare(si_id, seq, tscan_outdir)

    # Predict context plus scores
    contextplus_score_file = predict(mat_seq_file, ts_pred_file, utr_file)
    transcript_dict = get_transcript_dict(contextplus_score_file)
    if translation_dict is None:
        translation_dict = get_translation_dict(ref_seq_file)

    target_frame = parse(transcript_dict, translation_dict)

    out_fname = os.path.join(outdir, '%s.tsv' % args.id)
    target_frame.to_csv(out_fname, sep='\t', index=False)
    for f in [mat_seq_file, ts_pred_file, contextplus_score_file]:
        os.remove(f)
    return target_frame['GeneID', 'GeneName', 'GeneSynonyms'].drop_duplicates()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--id', help='siRNA identifier',
                        type=str)
    parser.add_argument('--seq',
                        help='siRNA antisense sequence w/o overhangs',
                        type=str)
    parser.add_argument('--tscan-outdir',
                        help='targetscan_60.pl output filepath',
                        type=str)
    parser.add_argument('--utr-file', help='UTR Seq filepath',
                        type=str)
    parser.add_argument('--ref-seq-file',
                        help='Reference file for transcript to gene lookup',
                        type=str)
    parser.add_argument('--outdir', help='Output directory filepath',
                        type=str)

    args = parser.parse_args()
    missing = (args.id is None or
               args.seq is None or
               args.tscan_outdir is None or
               args.utr_file is None or
               args.ref_seq_file is None or
               args.outdir is None)
    if missing:
        raise ValueError("Missing arguments")

    write_target_frame(args.id, args.seq, args.tscan_outdir,
                       args.utr_file, args.outdir, args.ref_seq_file)
