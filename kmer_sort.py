import random
import operator
import subprocess
from pathlib import Path
from argparse import ArgumentParser
from typing import Generator
import subprocess
from utils import *


def get_next_kmer(method: str, size: int, input: Path):

    if method == "random":
        return get_random_kmer(size=size, input=input)

    elif method == "occurence-max":
        return get_occ_kmer(method="max", input=input, size=size)

    elif method == "occurence-rate":
        return get_occ_kmer(method="rate", input=input, size=size)

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


@timer_func
def sort_by_kmer(input: Path, output: Path, size: int, method: str) -> None:
    """sort_by_kmer the read of the input file in the output file without the labels

    Args:
        input (str): input file
        output (str): output file
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

    with open_wipe_and_add(output) as file:
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


if __name__ == "__main__":
    pass
