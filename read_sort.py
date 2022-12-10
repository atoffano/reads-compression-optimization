from argparse import ArgumentParser
import tsne_sort, pca_sort, kmer_sort
from utils import *
import time


def argparser():
    """
    Parses arguments from command line
            Parameters:
                    Console arguments
                    -i / --input (str) : Input fasta file with headers
                    -m / --method (str): Method used to sort reads. Can be 'tsne_sort', 'pca_sort' or 'kmer_sort'
                    -d / --delete_output (bool) : Delete output file after compression ratio is computed
                    -c / --compare_to (str) : headerless fasta.gz file to compare to in order to compute compression ratio
                    -s / --size_kmer (int) : size of kmer used for kmer_sort
                    -cs / --chunk_size (int) : size of chunk used for tsne_sort and pca_sort
    """
    parser = ArgumentParser()

    # Add the arguments to the parser
    parser.add_argument("-i", "--input", required=False, help="")
    parser.add_argument("-d", "--delete_output", default=False, help="")
    parser.add_argument("-m", "--method", required=True, help="")
    parser.add_argument("-c", "--compare_to")
    parser.add_argument("-s", "--size_kmer", default=6, help="")
    parser.add_argument("-cs", "--chunk_size", default=1000000, help="")

    args = parser.parse_args()
    return args


def read_sort_main(
    input: Path,
    compare_to: Path,
    delete: bool = True,
    method: str = "kmer_sort",
    size: int = 6,
    chunk_size: int = 1000000,
):
    log = {}
    output = "_out.".join(input.split("."))

    t1 = time.time()

    if method == "kmer_sort":
        kmer_sort.sort_by_kmer(input=input, output=output, size=size, method="random")
    if method == "tsne_sort":
        tsne_sort.sort_by_tsne(input, output, int(chunk_size))
    if method == "pca_sort":
        pca_sort.sort_by_pca(input, output, int(chunk_size))

    log["time"] = time.time() - t1
    log["rate"] = monitor_gzip(output, compare_to)

    if delete:
        os.remove(output)
        os.remove(output + ".gz")

    return log


def main():
    args = argparser()
    args.delete_output = True if args.delete_output == "True" else False
    output = "_out.".join(args.input.split("."))

    if args.method == "kmer_sort":
        size = int(args.size_kmer)
        kmer_sort.sort_by_kmer(
            input=args.input, output=output, size=size, method="random"
        )

    if args.method == "tsne_sort":
        tsne_sort.sort_by_tsne(
            args.input, output, int(args.chunk_size)
        )

    if args.method == "pca_sort":
        pca_sort.sort_by_pca(
            args.input, output, int(args.chunk_size)
        )

    print(f"Compression ratio : {monitor_gzip(output, args.compare_to)}")

    if args.delete_output:
        os.remove(output)
        os.remove(output + ".gz")


if __name__ == "__main__":
    main()
    # Usage :
    # $ python read_sort.py -i data/ecoli_100Kb_reads_40x.fasta -m pca_sort -cs 40000 -c data/headerless/ecoli_100Kb_reads_40x.fasta.headerless.gz --delete_output True
    # python read_sort.py -i data/ecoli_100Kb_reads_120x.fasta -m kmer_sort -s 6 -c data/headerless/ecoli_100Kb_reads_120x.fasta.headerless.gz --delete_output True
