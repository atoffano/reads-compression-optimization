from argparse import ArgumentParser, Namespace
import tsne_sort, pca_sort, kmer_sort, chatGPT_sort
from utils import *
import time
from pathlib import Path


def argparser():
    parser = ArgumentParser()

    # Add the arguments to the parser
    parser.add_argument(
        "-i",
        "--infile",
        required=True,
        help="Path to source file to reorder",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=False,
        help="""Path to the output file if nothing provide output will """
        """be generated automatically in the same folder as source file""",
    )
    parser.add_argument(
        "-d",
        "--delete_output",
        default=False,
        help="if set as 'T' (True) the output file will be delete at the end",
    )
    parser.add_argument(
        "-c",
        "--compare_to",
        default=False,
        help="Path to to the file to compare the compressed output size with",
    )
    parser.add_argument(
        "-m",
        "--method",
        required=True,
        help="Method to use for reordering ('kmer_sort', 'pca_sort','tsne_sort',chatgpt_sort'",
    )
    parser.add_argument(
        "-s",
        "--size_kmer",
        default=6,
        help="Size of k-mer use in the 'kmer_sort' method",
    )
    parser.add_argument(
        "-cs",
        "--chunk_size",
        default=1000000,
        help="Size of chunk use in the 'pca_sort' and 'tsne_sort' methods",
    )
    parser.add_argument(
        "-cu",
        "--cutoff",
        default=0,
        help="Cutoff use in the 'kmer_sort' method if cutoff >= 0, no cutoff will be applied",
    )
    parser.add_argument(
        "-in",
        "--intervals_number",
        default=3,
        help="Numbers of intervals use in the 'kmer_sort' method",
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


def read_sort_main(
    infile: Path,
    compare_to: Path,
    delete: bool = True,
    method: str = "kmer_sort",
    size: int = 6,
    chunk_size: int = 1000000,
    intervals_number: int = 3,
    cutoff=0,
):
    log: dict = {}
    output = (
        Path("_organized.".join(infile.split(".")))
        if "." in infile
        else infile + "_organized"
    )
    infile = Path(infile)
    t1: float = time.time()

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
        os.remove(output.name + ".gz")

    return log


def main():
    """main function of the script

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
        chatGPT_sort.sort_by_minimizer(args.infile, args.output)

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
    # $ python read_sort.py -i data/ecoli_100Kb_reads_40x.fasta -m pca_sort -cs 40000 -c data/headerless/ecoli_100Kb_reads_40x.fasta.headerless.gz --delete_output True
    # python read_sort.py -i data/ecoli_100Kb_reads_120x.fasta -m kmer_sort -s 6 -c data/headerless/ecoli_100Kb_reads_120x.fasta.headerless.gz --delete_output True
