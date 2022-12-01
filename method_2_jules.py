import random
import operator
import subprocess
from pathlib import Path
from argparse import ArgumentParser
from typing import Generator
import subprocess
from utils import fasta_reader, get_sequence, gzip_out, open_wipe_and_add
from profiling import monitor, monitor_gzip


def get_next_kmer(method: str, size: int, input: Path):

    if method == "random":
        return get_random_kmer(size=size, input=input)

    elif method == "occurence-max":
        return get_occ_kmer(method="max")

    elif method == "occurence-rate":
        return get_occ_kmer(method="rate")

    elif method == "first-found":
        return get_first_kmer(input=input, size=size)

    else:
        raise "Method do not exist"


def get_first_kmer(input, size: int):
    return get_sequence(next(fasta_reader(input)))[0:size]


def get_random_kmer(size: str, input: Path) -> str:

    seq: str = get_sequence(next(fasta_reader(input)))
    start: int = random.randint(0, len(seq) - size)
    return seq[start : start + size]


def get_occ_kmer(method: str, input: Path, size: int) -> str:

    reader: Generator = fasta_reader(filename=input)
    kmers_score: dict = {}
    for read in reader:
        kmer_found_in_seq = False
        seq = get_sequence(read)
        for i in range(len(seq) - size):
            kmer = seq[i : i + size]

            if not kmer in kmers_score:
                kmers_score[kmer] = [1, 1]
                kmer_found_in_seq = True
            else:
                kmers_score[kmer][0] += 1
                # kmers_score[kmer] += 1
                if not kmer_found_in_seq:
                    kmers_score[kmer][1] += 1

    if method == "rate":
        for key in kmers_score:
            kmers_score[key] = kmers_score[key][0] / kmers_score[key][1]

    elif method == "max":
        for key in kmers_score:
            kmers_score[key] = kmers_score[key][0]

    return max(kmers_score.items(), key=operator.itemgetter(1))[0]


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


def get_kmer_dict(kmer: str, input: Path, size: int, method: str) -> dict:

    kmer_dict: dict = {}
    round_limit: int = 100
    breaking: bool = False

    # a dict is build for the first kmer and then we pick another kmer to sort the
    # reads in which the previous kmer had 0 occurences
    for _ in range(round_limit):

        # build the occurence dict for the current kmer
        kmer_dict[kmer]: dict = get_kmer_count_dict(kmer, input)

        # if all the reads has been sorted we end the loop
        if not "0" in kmer_dict[kmer]:
            breaking: bool = True
            break

        else:
            with open_wipe_and_add("method_2_temp") as file:
                for read in kmer_dict[kmer]["0"]:
                    file.write(read + "\n")

            last_zero: list = kmer_dict[kmer]["0"]
            del kmer_dict[kmer]["0"]

        input: Path = "method_2_temp"
        kmer: str = get_next_kmer(method=method, size=size, input=input)

    if not breaking:
        # if we reach 100 kmer we keep the rest of the reads as an unsorted list
        kmer_dict["LAST"]: dict = {"0": last_zero}
    return kmer_dict


def reorder(input: Path, out_file: Path, size: int, method: str) -> None:
    """reorder the read of the input file in the output file without the labels

    Args:
        input (str): input file
        out_file (str): output file
        kmer (str, optional): first kmer to use may be random. Defaults to "".
        size (str, optional): if first kmer is random set size of kmer to use. Defaults to "".
    """
    first_kmer: str = get_next_kmer(method=method, size=size, input=input)

    kmer_dict: dict = get_kmer_dict(
        method=method, kmer=first_kmer, input=input, size=size
    )

    kmer_list: list = [kmer for kmer in kmer_dict.keys()]
    kmer_list.sort()
    sorted_kmer_dict: dict = {}

    for kmer in kmer_list:
        sorted_kmer_dict[kmer] = kmer_dict[kmer]

    with open_wipe_and_add(out_file) as file:
        for kmer in sorted_kmer_dict:
            for count in sorted_kmer_dict[kmer]:
                for read in sorted_kmer_dict[kmer][count]:
                    file.write(get_sequence(read) + "\n")


def erro_handling(input: Path, size: int) -> None:
    """Handle error regarding the kmer or size given in arguments

    Args:
        input (Path): input file
        kmer (str): given kmer
        size (int): given size

    Raises:
        ValueError: kmer size check
        ValueError: positive and non-null size check
        ValueError: size limit check
    """
    reader: Generator = fasta_reader(filename=input)
    read = next(reader)

    if size <= 0:
        raise ValueError("Size cannot be equal or lesser than 0")

    if size > len(get_sequence(read)):
        raise ValueError("Size cannot be greater than reads size")


def launch_method_2(
    size: int,
    method: str,
    input: Path,
    output: Path,
    delete_output: bool = False,
) -> dict:
    """function for launching the method inside a python script

    Args:
        input (Path): input file path
        output (Path): output file path (may be deleted)
        delete_output (bool, optional): if True delete output file after rate computation. Defaults to False.
        kmer (str, optional): the first kmer to use. Defaults to "random".
        size (int, optional): if kmer is empty set the kmer size to use. Defaults to 6.

    Returns:
        dict: logs of the operation
    """

    random.seed(42)
    erro_handling(input=input, size=size)
    # Reorder and get monitoring values
    monitoring_values = monitor(
        reorder(input=input, out_file=output, size=size, method=method)
    )

    # Compress the files
    in_gz = "data/ori/" + input.split("data/")[1] + ".gz"
    out_gz = gzip_out(output)

    # get compression rate
    rate = monitor_gzip(out_gz, in_gz)

    # Delete the file if -d is on
    if delete_output:
        bashCommand = f"rm -rf {out_gz}"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        out, error = process.communicate()

    monitoring_values["file"] = input
    monitoring_values["kmer-size"] = size
    monitoring_values["rate"] = rate

    return monitoring_values


def main() -> None:
    """main function to launch the method with command line"""

    # argparse
    parser = ArgumentParser(
        prog="method_2_jules.py",
        description="",
        epilog="Text at the bottom of help",
    )

    parser.add_argument("-i", "--input")
    parser.add_argument("-d", "--delete_output", default=False)
    parser.add_argument("-o", "--output")
    parser.add_argument("-s", "--size_kmer")
    parser.add_argument("-m", "--method")

    args = parser.parse_args()
    size = int(args.size_kmer)
    delete_output = False if args.delete_output == "False" else True

    # Check for errors
    erro_handling(input=args.input, size=size)

    # Reorder and get monitoring values
    print("reordering...\n")
    monitoring_values = monitor(
        reorder(
            input=args.input,
            out_file=args.output,
            size=size,
            method=args.method,
        )
    )

    # Print monitoring values
    for key in monitoring_values:
        print(f"{key} : {monitoring_values[key]}")

    # Compress the files
    in_gz = "data/ori/" + args.input.split("data/")[1] + ".gz"
    out_gz = gzip_out(args.output)

    # Print compression rate
    compression_rate = monitor_gzip(out_gz, in_gz)

    # Delete the file if -d is on
    if delete_output:
        bashCommand = f"rm -rf {out_gz}"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        out, error = process.communicate()

    # write log in file
    with open_wipe_and_add("log_method2.txt") as file:
        for key in monitoring_values:
            file.write(f"{key} : {monitoring_values[key]}\n")

        file.write(f"compression rate : {compression_rate}")


if __name__ == "__main__":

    main()
