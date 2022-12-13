from argparse import ArgumentParser
import tsne_sort, pca_sort, kmer_sort, chatGPT_sort
from utils import *
import time


def argparser():
    """
    Parses arguments from command line
            Parameters:
                    Console arguments
                    -i / --infile (str) : infile fasta file with headers
                    -m / --method (str): Method used to sort reads. Can be 'tsne_sort', 'pca_sort' or 'kmer_sort'
                    -d / --delete_output (bool) : Delete output file after compression ratio is computed
                    -c / --compare_to (str) : headerless fasta.gz file to compare to in order to compute compression ratio
                    -s / --size_kmer (int) : size of kmer used for kmer_sort
                    -cs / --chunk_size (int) : size of chunk used for tsne_sort and pca_sort
    """
    parser = ArgumentParser()

    # Add the arguments to the parser
    parser.add_argument("-i", "--infile", required=False, help="")
    parser.add_argument("-d", "--delete_output", default="True", help="")
    parser.add_argument("-m", "--method", required=True, help="")
    parser.add_argument("-c", "--compare_to")
    parser.add_argument("-s", "--size_kmer", default=6, help="")
    parser.add_argument("-cs", "--chunk_size", default=0, help="")
    parser.add_argument("-cu", "--cutoff", default=0, help="")
    parser.add_argument("-in", "--intervals_number", default=3, help="")

    args = parser.parse_args()
    return args


def read_sort_main(
    infile: Path,
    compare_to: Path,
    delete: bool = True,
    method: str = "kmer_sort",
    size: int = 6,
    chunk_size: int = 0,
    intervals_number: int = 3,
    cutoff=0,
):
    log = {}
    output = "_out.".join(infile.split("."))

    t1 = time.time()

    if method == "kmer_sort":
        kmer_sort.sort_by_kmer(
            infile=infile,
            output=output,
            size=size,
            intervals_number=intervals_number,
            cutoff=cutoff,
        )
    elif method == "tsne_sort":
        tsne_sort.sort_by_tsne(infile, output, int(chunk_size))
    elif method == "pca_sort":
        pca_sort.sort_by_pca(infile, output, int(chunk_size))

    log["time"] = time.time() - t1
    log["rate"] = monitor_gzip(output, compare_to)

    if delete:
        os.remove(output)
        os.remove(output + ".gz")

    return log


def main():
    args = argparser()
    args.delete_output = True if args.delete_output == "True" else False
    output = "_out.".join(args.infile.split("."))

    if args.method == "kmer_sort":
        size = int(args.size_kmer)
        cutoff = int(args.cutoff)
        intervals_number = int(args.intervals_number)
        kmer_sort.sort_by_kmer(
            infile=args.infile,
            output=output,
            size=size,
            intervals_number=intervals_number,
            cutoff=cutoff,
        )

    if args.method == "tsne_sort":
        tsne_sort.sort_by_tsne(args.infile, output, int(args.chunk_size))

    if args.method == "pca_sort":
        pca_sort.sort_by_pca(args.infile, output, int(args.chunk_size))

    if args.method == "chatgpt_sort":
        chatGPT_sort.sort_by_minimizer(args.infile, output)

    print(f"Compression ratio : {monitor_gzip(output, args.compare_to)}")

    if str(args.delete_output) == "True":
        os.remove(output)
        os.remove(output + ".gz")


if __name__ == "__main__":
    main()
    # Usage :
    # $ python read_sort.py -i data/ecoli_100Kb_reads_40x.fasta -m pca_sort -cs 40000 -c data/headerless/ecoli_100Kb_reads_40x.fasta.headerless.gz --delete_output True
    # python read_sort.py -i data/ecoli_100Kb_reads_120x.fasta -m kmer_sort -s 6 -c data/headerless/ecoli_100Kb_reads_120x.fasta.headerless.gz --delete_output True
