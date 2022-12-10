import random
from copy import deepcopy

# import operator
from pathlib import Path
from typing import Generator
from utils import *


# def get_next_kmer(method: str, size: int, infile: Path) -> str:

#     if method == "random":
#         return get_random_kmer(size=size, infile=infile)

#     elif method == "occurence-max":
#         return get_occ_kmer(method="max", infile=infile, size=size)

#     elif method == "occurence-rate":
#         return get_occ_kmer(method="rate", infile=infile, size=size)

#     elif method == "first-found":
#         return get_first_kmer(infile=infile, size=size)

#     else:
#         raise "Method do not exist"


# def get_first_kmer(infile: str, size: int) -> str:
#     return get_sequence(next(fasta_reader(infile)))[0:size]


# def get_occ_kmer(method: str, infile: Path, size: int) -> str:

#     reader: Generator = fasta_reader(filename=infile)
#     kmers_score: dict = {}
#     for read in reader:
#         kmer_found_in_seq = False
#         seq = get_sequence(read)
#         for i in range(len(seq) - size):
#             kmer = seq[i : i + size]

#             if not kmer in kmers_score:
#                 kmers_score[kmer] = [1, 1]
#                 kmer_found_in_seq = True
#             else:
#                 kmers_score[kmer][0] += 1
#                 # kmers_score[kmer] += 1
#                 if not kmer_found_in_seq:
#                     kmers_score[kmer][1] += 1

#     if method == "rate":
#         for key in kmers_score:
#             kmers_score[key] = kmers_score[key][0] / kmers_score[key][1]

#     elif method == "max":
#         for key in kmers_score:
#             kmers_score[key] = kmers_score[key][0]

#     return max(kmers_score.items(), key=operator.itemgetter(1))[0]
# def split_in_chunks(infile: Path, split_size: int) -> list:

#     line_no = len([1 for line in fasta_reader(infile) if line]) * 2
#     if line_no > 240000 or split_size == 0:
#         chunk_files = []
#         sep = line_no // 240000
#         chunk_sizes = [int(line_no / sep) for _ in range(sep)]
#         if line_no % 240000 != 0:
#             chunk_sizes.append(line_no % 240000)

#         reader = fasta_reader(infile)
#         for i, chunksize in enumerate(chunk_sizes):
#             chunk_name = f"tempfile/kmer_sort_c{i}_temp"
#             chunk_files.append(chunk_name)
#             with open_wipe_and_add(chunk_name) as chunk:
#                 for _ in range(chunksize):
#                     try:
#                         chunk.write(next(reader) + "\n")
#                     except StopIteration:
#                         pass

#         return chunk_files

#     else:
#         return [infile]


def get_random_kmer(size: str, infile: Path) -> str:

    seq: str = get_sequence(next(fasta_reader(infile)))
    start: int = random.randint(0, len(seq) - size)
    return seq[start : start + size]


def get_identifier(read: str, kmer: str, intervals: list) -> str:
    seq = get_sequence(read)
    kmer_positions = ""
    last_kmer_index = 0
    move_index = 1
    while move_index != 0:  # if the result of the find() is -1 we stop the loop

        move_index = seq[last_kmer_index:].find(kmer) + 1

        # a = input()
        if move_index > 0:
            last_kmer_index += move_index

            if last_kmer_index >= intervals[-1]:
                kmer_positions += str(len(intervals))
            else:
                kmer_positions += str(
                    [j for j, int in enumerate(intervals) if last_kmer_index > int][-1]
                    + 1
                )

    return str(len(kmer_positions)) + kmer_positions


def get_kmer_count_dict(
    kmer: str, current_file: Path, next_file: Path, intervals_number: int, seqlen: int
) -> dict:
    """Build the dict with the num of occurence of the given kmer as key and the lists
       of the correspondings reads as value

    Args:
        kmer (str): the given kmer
        file (Path): the infile file path

    Returns:
        dict: dict with num of occurence as keys and reads lists as values
    """
    not_sorted: list = []
    kmer_count_dict: dict = {}
    reader: Generator = fasta_reader(current_file)
    intervals = [i * (seqlen / intervals_number) for i in range(intervals_number)]
    with open_wipe_and_add(next_file) as file:
        for read in reader:

            identifier: str = get_identifier(read=read, kmer=kmer, intervals=intervals)
            if identifier == "0":
                not_sorted.append(read)
                file.write(read + "\n")
            else:
                if not identifier in kmer_count_dict:
                    kmer_count_dict[identifier] = [read]
                else:
                    kmer_count_dict[identifier].append(read)

    return kmer_count_dict, not_sorted


def get_kmer_dict(
    kmer: str,
    infile: Path,
    original_size: int,
    intervals_number: int,
    seqlen: int,
    cutoff: int,
) -> dict:

    kmer_dict: dict = {}
    breaking: bool = False
    size = original_size
    current_file = infile
    next_file = "tempfile/kmer_sort_no_init"
    # a dict is build for the first kmer and then we pick another kmer to sort the
    # reads in which the previous kmer had 0 occurences

    for round in range(cutoff):

        # build the occurence dict for the current kmer
        kmer_dict[kmer], not_sorted = get_kmer_count_dict(
            kmer=kmer,
            current_file=current_file,
            next_file=next_file,
            intervals_number=intervals_number,
            seqlen=seqlen,
        )
        if current_file != infile:
            os.remove(current_file)

        current_file = next_file
        next_file = current_file.split("no")[0] + f"no{str(round)}"

        # if all the reads has been sorted we end the loop
        if not_sorted == []:
            breaking: bool = True
            break

        else:
            last_zero: list = deepcopy(not_sorted)

        kmer: str = get_random_kmer(size=size, infile=current_file)
        current_size = size
        while kmer in kmer_dict:
            kmer: str = get_random_kmer(size=size, infile=current_file)
            size += 1
        size = current_size

    if not breaking:
        # if we reach 1000 kmer we keep the rest of the reads as an unsorted list
        print("reaching 1000")
        kmer_dict["LAST"]: dict = {"0": last_zero}

    return kmer_dict


@timer_func
def sort_by_kmer(
    infile: Path, outfile: Path, size: int, intervals_number: int, cutoff: int
) -> None:
    """sort_by_kmer the read of the infile file in the outfile file without the labels

    Args:
        infile (str): infile file
        outfile (str): outfile file
        size (str, optional): if first kmer is random set size of kmer to use. Defaults to "".
    """
    random.seed(42)
    seqlen = len(get_sequence(next(fasta_reader(infile))))
    first_kmer: str = get_random_kmer(size=size, infile=infile)
    kmer_dict: dict = get_kmer_dict(
        kmer=first_kmer,
        infile=infile,
        original_size=size,
        intervals_number=intervals_number,
        seqlen=seqlen,
        cutoff=cutoff,
    )

    with open_wipe_and_add(outfile) as file:
        for kmer in kmer_dict:
            for count in kmer_dict[kmer]:
                for read in kmer_dict[kmer][count]:
                    file.write(get_sequence(read) + "\n")


def erro_handling(infile: Path, size: int) -> None:
    """Handle error regarding the kmer or size given in arguments

    Args:
        infile (Path): infile file
        kmer (str): given kmer
        size (int): given size

    Raises:
        ValueError: kmer size check
        ValueError: positive and non-null size check
        ValueError: size limit check
    """
    reader: Generator = fasta_reader(filename=infile)
    read = next(reader)

    if size <= 0:
        raise ValueError("Size cannot be equal or lesser than 0")

    if size > len(get_sequence(read)):
        raise ValueError("Size cannot be greater than reads size")


if __name__ == "__main__":
    pass

# python read_sort.py -i data/ecoli_100Kb_reads_120x.fasta -m kmer_sort -s 8 -cu 3 -c data/headerless/ecoli_100Kb_reads_120x.fasta.headerless.gz -d True
# python read_sort.py -i data/humch1_1Mb_reads_120x.fasta -m kmer_sort -s 8 -cu 3 -c data/headerless/humch1_1Mb_reads_120x_headerless.fasta.gz -d True
#
# 2 : 1.13
# 3 : 1.16
# 4 : 1.23
#
# 6 : 1.78
# 6 cs 240_000 : 1.72
# 6 cs 240_000 + ps : 1.72
# 6 cs 300_000 : 1.72
# 7 : 1.62
