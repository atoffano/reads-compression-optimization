from argparse import ArgumentParser, Namespace
import tsne_sort, pca_sort, kmer_sort, rollinghash_sort
from utils import *
import time
from pathlib import Path


def argparser():
    """Parse the arguments passed to the script"""

    # Create the parser
    parser = ArgumentParser()

    # Add the arguments to the parser
    parser.add_argument("-i", "--infile", required=True, help="Input file path")
    parser.add_argument(
        "-o",
        "--output",
        required=False,
        help="Output file path, if not provided output will be save in the script file folder",
    )
    parser.add_argument(
        "-d",
        "--delete_output",
        default=False,
        help="""passed 'T' (True) to delete any output after the algorithm"""
        """. This only should be used to monitor performance""",
    )
    parser.add_argument(
        "-m",
        "--method",
        required=True,
        help="Method to use between 'tsne_sort', 'pca_sort','kmer_sort','rollinghash_sort'",
    )
    parser.add_argument(
        "-c",
        "--compare_to",
        default=False,
        help="File path to compare the compressed output file with",
    )
    parser.add_argument(
        "-s", "--size_kmer", default=6, help="Kmer size to use in kmer_sort method"
    )
    parser.add_argument(
        "-cs",
        "--chunk_size",
        default=0,
        help="""Chunk size to use in pca_sort and tsne_sort method, by default chunk size"""
        """will be equal to number of read in the file""",
    )
    parser.add_argument(
        "-cu",
        "--cutoff",
        default=0,
        help="Cutoff to use in kmer_sort method, by default there is no cutoff",
    )
    parser.add_argument(
        "-in",
        "--intervals_number",
        default=3,
        help="Number of intervals to use in the kmer_sort method",
    )

    args = parser.parse_args()

    return argparse_parser(args)


def argparse_parser(args: Namespace) -> Namespace:
    """Reparse the argparse argument

    Args:
        args (Namespace): argparse arguments

    Raises:
        ValueError: raise error if -d argument is different than 'T'

    Returns:
        Namespace: the parsed Namespace instance
    """
    args.infile = Path(args.infile)
    args.compare_to = Path(args.compare_to) if args.compare_to else False
    if not args.output:
        args.output = (
            Path("_organized.".join(args.infile.name.split(".")))
            if args.infile.name.count(".") == 1
            else args.infile.name + "_organized"
        )
    else:
        args.output = Path(args.output)

    if args.delete_output == "T":
        args.delete_output = True
    elif args.delete_output != False:
        raise ValueError(
            "-d / --delete_output arg value must be 'T' (True) or 'F' (False)"
        )
    args.size_kmer = int(args.size_kmer)
    args.chunk_size = int(args.chunk_size)
    args.cutoff = int(args.cutoff)
    args.intervals_number = int(args.intervals_number)

    return args

def main():
    """Main function of the script. Calls the argument parser, then call the sorting algorithm.

    Raises:
        ValueError: raise error if the method provided does not exist
    """
    args = argparser()

    if args.method == "kmer_sort":
        kmer_sort.sort_by_kmer(
            infile=args.infile,
            output=args.output,
            size=args.size_kmer,
            intervals_number=args.intervals_number,
            cutoff=args.cutoff,
        )

    elif args.method == "tsne_sort":
        tsne_sort.sort_by_tsne(args.infile, args.output, args.chunk_size)

    elif args.method == "pca_sort":
        pca_sort.sort_by_pca(args.infile, args.output, args.chunk_size)

    elif args.method == "chatgpt_sort":
        rollinghash_sort.sort_by_minimizer(args.infile, args.output)

    else:
        raise ValueError(
            "Method does not exist try 'kmer_sort','pca_sort','tsne_sort', or 'chatgpt_sort'"
        )
    if args.compare_to:
        print(f"Compression ratio : {monitor_gzip(args.output, args.compare_to)}")

    if (
        args.delete_output
    ):  # this should only be use to test the performance of the sorting algorithm
        os.remove(args.output)
        if args.compare_to:
            os.remove(args.output.name + ".gz")


if __name__ == "__main__":
    main()
    # Usage :
    # python read_sort.py -i data/ecoli_100Kb_reads_5x.fasta -m kmer_sort -d T
    # python read_sort.py -i data/ecoli_100Kb_reads_5x.fasta -m pca_sort -d T
    # python read_sort.py -i data/ecoli_100Kb_reads_5x.fasta -m tsne_sort -d T
    # python read_sort.py -i data/ecoli_100Kb_reads_5x.fasta -m chatgpt_sort -d T
