from argparse import ArgumentParser
import tsne_sort, pca_sort, kmer_sort
from utils import *

def argparser():
    '''
    Parses arguments from command line
            Parameters:
                    Console arguments
                    -i / --input (str) : Input fasta file with headers
                    -m / --method (str): Method used to sort reads. Can be 'tsne_sort', 'pca_sort' or 'kmer_sort'
                    -d / --delete_output (bool) : Delete output file after compression ratio is computed
                    -c / --compare_to (str) : headerless fasta.gz file to compare to in order to compute compression ratio
                    -s / --size_kmer (int) : size of kmer used for kmer_sort
                    -cs / --chunk_size (int) : size of chunk used for tsne_sort and pca_sort
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
    # Usage :
    # $ python read_sort.py -i data/ecoli_100Kb_reads_40x.fasta -m pca_sort -cs 40000 -c data/headerless/ecoli_100Kb_reads_40x.fasta.headerless.gz --delete_output True