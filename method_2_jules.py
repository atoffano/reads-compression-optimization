import random
from pathlib import Path
from argparse import ArgumentParser
from typing import Generator

from utils import fasta_reader, get_sequence

# géré le cas ou l'utilisateur donne un kmer trop long, un size négatif, un size trop long
# trouvé une méthode pour que la sélection des kmers soit smart et pas aléatoire


def get_first_kmer(kmer: str, size: str, in_file: Path) -> str:

    if kmer != None:
        return kmer  # if a kmer is give as parameter use it as first kmer
    else:
        try:
            size = int(size)
        except ValueError:
            raise ValueError("kmer size (-s, --size-mer) must be an integer")

        # if no kmer given in args, select a random kmer in the first read
        reader: Generator = fasta_reader(in_file)
        read_seq: str = get_sequence(next(reader))
        kmer_start: int = random.randint(0, len(read_seq) - size)

        return read_seq[kmer_start : kmer_start + size]


def get_kmer_count_dict(kmer: str, file: Path) -> dict:
    """Build the dict with the num of occurence of the given kmer as key and the lists
       of the correspondings reads as value

    Args:
        kmer (str): the given kmer
        file (Path): the input file path

    Returns:
        dict: dict with num of occurence as keys and reads lists as values
    """
    kmer_count_dict: dict = {}
    reader: Generator = fasta_reader(file)
    for read in reader:
        kmer_count: str = str(get_sequence(read).count(kmer))
        if not kmer_count in kmer_count_dict:
            kmer_count_dict[kmer_count] = [read]
        else:
            kmer_count_dict[kmer_count].append(read)

    return kmer_count_dict


def get_kmer_dict(kmer: str, in_file: Path) -> dict:

    kmer_dict: dict = {}
    round_limit: int = 100
    breaking: bool = False

    # a dict is build for the first kmer and then we pick another kmer to sort the
    # reads in which the previous kmer had 0 occurences
    for _ in range(round_limit):

        # build the occurence dict for the current kmer
        kmer_dict[kmer]: dict = get_kmer_count_dict(kmer, in_file)

        # if all the reads has been sorted we end the loop
        if not "0" in kmer_dict[kmer]:
            breaking: bool = True
            break

        else:
            # pick a random kmer in the first sequence of the reads that had 0
            # occurences of the previous kmer
            new_kmer_read: str = get_sequence(kmer_dict[kmer]["0"][0])
            index: int = random.randint(0, len(new_kmer_read) - len(kmer))
            new_kmer: str = new_kmer_read[index : index + len(kmer)]

            # write the reads with 0 occurence of the previous kmer as new input file
            with open("method_2_temp", "w") as wipe:
                wipe.write("")
            with open("method_2_temp", "a") as file:
                for read in kmer_dict[kmer]["0"]:
                    file.write(read + "\n")

            last_zero: list = kmer_dict[kmer]["0"]
            del kmer_dict[kmer]["0"]

        in_file: Path = "method_2_temp"
        kmer: str = new_kmer

    if not breaking:
        # if we reach 100 kmer we keep the rest of the reads as an unsorted list
        kmer_dict["LAST"]: dict = {"0": last_zero}
    return kmer_dict


def write_outfile(in_file: str, out_file: str, kmer: str = "", size: str = "") -> None:
    in_file = Path(in_file)
    out_file = Path(out_file)

    first_kmer: str = get_first_kmer(kmer=kmer, size=size, in_file=in_file)
    kmer_dict: dict = get_kmer_dict(first_kmer, in_file)
    kmer_list = [kmer for kmer in kmer_dict.keys()]
    kmer_list.sort()
    sorted_kmer_dict = {}

    for kmer in kmer_list:
        sorted_kmer_dict[kmer] = kmer_dict[kmer]

    with open(out_file, "w") as wipe:
        wipe.write("")
    with open(out_file, "a") as file:
        for kmer in sorted_kmer_dict:
            for count in sorted_kmer_dict[kmer]:
                for read in sorted_kmer_dict[kmer][count]:
                    file.write(read + "\n")


if __name__ == "__main__":
    parser = ArgumentParser(
        prog="method_2_jules.py",
        description="",
        epilog="Text at the bottom of help",
    )

    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    parser.add_argument("-k", "--kmer")
    parser.add_argument("-s", "--size_kmer")

    args = parser.parse_args()

    if args.kmer == None and args.size_kmer == None:
        args.size_kmer == "6"

    arg_log_kmer: str = "random" if args.kmer == None else args.kmer
    arg_log_size: int = args.size_kmer if args.kmer == None else len(args.kmer)

    print(
        f"Parameters are :\n\tFirst kmer : {arg_log_kmer}\n\tKmer size : {arg_log_size}"
    )

    write_outfile(
        in_file=args.input, out_file=args.output, kmer=args.kmer, size=args.size_kmer
    )
