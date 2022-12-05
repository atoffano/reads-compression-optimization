from argparse import ArgumentParser
import tsne_sort, pca_sort, kmer_sort
from utils import *

def argparser():
    '''
    Parses arguments from command line
            Parameters:
                    Console arguments
                    -i / --input (str) : Input path containing files to convert. Handles folder containing multiple datasets.
                    -m / --method (str): Output path for both standardized and converted formats.
                    -s / --score (str) : Input format.
                    -e / --evalset (str) : Output format.
                    -o / --original (str) : Output format.
    '''
    parser = ArgumentParser()

    # Add the arguments to the parser
    parser.add_argument("-i", "--input", required=False,
    help="")
    parser.add_argument("-d", "--delete_output", default=False,
    help="")
    parser.add_argument("-m", "--method", required=True,
    help="")
    parser.add_argument("-c", "--compare_to")
    parser.add_argument("-s", "--size_kmer",
    help="")
    parser.add_argument("-cs", "--chunk_size", default=1000000,
    help="")

    args = parser.parse_args()
    return args

def main():
    args = argparser()

    if args.method == 'kmer_sort':
        size = int(args.size_kmer)

    if args.method == 'tsne_sort':
        tsne_sort.sort_by_tsne(args.input, "out_x.fasta", int(args.chunk_size))

    if args.method == 'pca_sort':
        pca_sort.sort_by_pca(args.input, "out_x.fasta", int(args.chunk_size))
    
    print(f'Compression ratio : {monitor_gzip("out_x.fasta", args.compare_to)}')

    if args.delete_output:
        os.remove('out_x.fasta')
        os.remove('out_x.fasta.gz')

if __name__ == "__main__":
    main()
    # Utilization example : $ python read_sort.py -i data/ecoli_100Kb_reads_40x.fasta -m pca_sort -cs 40000 -c data/headerless/ecoli_100Kb_reads_40x.fasta.headerless.gz --delete_output True