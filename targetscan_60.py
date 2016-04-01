import os
import subprocess
import multiprocessing as mp


def run_targetscan(seed, utr_seq_file, seeds_dir, out_dir):
    outfile = os.path.join(out_dir, 'tscan.%s.tsv' % seed)
    if not os.path.isfile(outfile):
        script = """
                TSCAN=$(which targetscan_60.pl)
                perl $TSCAN {0}/{1}.tsv {2} {3}
                """.format(seeds_dir, seed, utr_seq_file, outfile)
        proc = subprocess.Popen(['bash', '-c', script],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                stdin=subprocess.PIPE)
        stdout, stderr = proc.communicate()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--utr-seq-file',
                        help="Filepath to 3' UTR sequence file",
                        type=str)
    parser.add_argument('--seeds-directory',
                        help="Filepath to seed sequence directory",
                        type=str)
    parser.add_argument('--out-directory', help="Filepath to output directory",
                        type=str)
    parser.add_argument('--n-workers',
                        help='Number of parallel processes to use',
                        type=int, default=12)

    args = parser.parse_args()
    seed_sequences = map(lambda x: os.path.splitext(x)[0],
                         os.listdir(args.seeds_directory))

    def foo(seed):
        run_targetscan(seed, args.utr_seq_file,
                       args.seeds_directory,
                       args.out_directory)
    if args.n_workers == 1:
        for seed in seed_sequences:
            foo(seed)
    else:
        p = mp.Pool(args.n_workers)
        p.map(foo, seed_sequences)
        p.terminate()
