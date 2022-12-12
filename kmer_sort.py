import random
from typing import Generator, Tuple
from utils import *
from pathlib import Path


def get_random_kmer(size: int, infile: Path) -> str:
    """return a random k-mer picked in the first sequence of the current source file

    Args:
        size (int): size of the kmer to generate
        infile (Path): source file path

    Returns:
        str: _description_
    """
    seq: str = get_sequence(next(fasta_reader(infile)))
    start: int = random.randint(0, len(seq) - size)
    return seq[start : start + size]


def get_identifier(read: str, kmer: str, intervals: list) -> str:
    seq: str = get_sequence(read)
    kmer_positions: str = ""
    last_kmer_index: int = 0
    move_index: int = 1
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


def get_identifiers_dict(
    kmer: str,
    current_file: str,
    next_file: str,
    intervals_number: int,
    seqlen: int,
) -> Tuple[dict, list]:
    """Build the dict with the num of occurence of the given kmer as key and the lists
       of the correspondings reads as value

    Args:
        kmer (str): the given kmer
        file (str): the infile file str

    Returns:
        dict: dict with num of occurence as keys and reads lists as values
    """
    not_sorted: list = []
    identifiers_dict: dict = {}
    reader: Generator = fasta_reader(current_file)
    intervals: list = [i * (seqlen / intervals_number) for i in range(intervals_number)]
    with open_wipe_and_add(next_file) as file:
        for read in reader:

            identifier: str = get_identifier(read=read, kmer=kmer, intervals=intervals)
            if identifier == "0":
                not_sorted.append(get_sequence(read))
                file.write(read + "\n")
            else:
                if not identifier in identifiers_dict:
                    identifiers_dict[identifier] = [get_sequence(read)]
                else:
                    identifiers_dict[identifier].append(get_sequence(read))

    return identifiers_dict, not_sorted


def sorting(
    kmer: str,
    infile: str,
    output: str,
    original_size: int,
    intervals_number: int,
    seqlen: int,
    cutoff: int,
) -> dict:

    used_kmer: list = []
    size: int = original_size
    current_file: str = infile
    next_file: str = "tempfile/kmer_sort_no_init"
    with open_wipe_and_add(output) as output_file:
        for round in range(cutoff):

            used_kmer.append(kmer)

            identifiers_dict, not_sorted = get_identifiers_dict(
                kmer=kmer,
                current_file=current_file,
                next_file=next_file,
                intervals_number=intervals_number,
                seqlen=seqlen,
            )

            for identifier in identifiers_dict:
                for seq in identifiers_dict[identifier]:
                    output_file.write(seq + "\n")

            # if all the reads has been sorted we end the loop
            if not not_sorted:
                break

            if current_file != infile:
                os.remove(current_file)

            current_file = next_file
            next_file = current_file.split("no")[0] + f"no{str(round)}"

            current_size = size
            while kmer in used_kmer:
                kmer: str = get_random_kmer(size=size, infile=current_file)
                size += 1
            size = current_size

        if not_sorted:
            # if we reach the cutoff we keep the rest of the reads as an unsorted list
            print(f"reaching {cutoff}")
            for seq in not_sorted:
                output_file.write(seq + "\n")

    if current_file != infile:
        os.remove(current_file)
        os.remove(next_file)


@timer_func
def sort_by_kmer(
    infile: str, output: str, size: int, intervals_number: int, cutoff: int
) -> None:
    """sort_by_kmer the read of the infile file in the output file without the labels

    Args:
        infile (str): infile file
        output (str): output file
        size (str, optional): if first kmer is random set size of kmer to use. Defaults to "".
    """
    random.seed(42)

    erro_handling(
        infile=infile,
        output=output,
        size=size,
        intervals_number=intervals_number,
        cutoff=cutoff,
    )
    cutoff = len(list(fasta_reader(infile))) if cutoff <= 0 else cutoff
    seqlen = len(get_sequence(next(fasta_reader(infile))))
    first_kmer: str = get_random_kmer(size=size, infile=infile)
    sorting(
        kmer=first_kmer,
        infile=infile,
        output=output,
        original_size=size,
        intervals_number=intervals_number,
        seqlen=seqlen,
        cutoff=cutoff,
    )


def erro_handling(
    infile: str, output: str, size: int, intervals_number: int, cutoff: int
) -> None:
    """Handle error regarding the kmer or size given in arguments

    Args:
        infile (str): infile file
        kmer (str): given kmer
        size (int): given size

    Raises:
        ValueError: kmer size check
        ValueError: positive and non-null size check
        ValueError: size limit check
    """
    reader: Generator = fasta_reader(filename=infile)
    read = next(reader)

    if size > len(get_sequence(read)):
        raise ValueError("Size cannot be greater than reads size")
    if size <= 1:
        raise ValueError("Size must be greater than 1")
    if intervals_number > len(get_sequence(read)):
        raise ValueError("interval number cannot be greater than reads size")
    if intervals_number < 1:
        raise ValueError("Intervals number must greater than 0")
    int_size = len(get_sequence(read)) / intervals_number
    if int_size <= 10:
        print(
            f"Warning : interval size are close to original kmer_size it could reduce sorting performance"
        )


if __name__ == "__main__":
    pass

# python read_sort.py -i data/ecoli_100Kb_reads_120x.fasta -m kmer_sort -c data/headerless/ecoli_100Kb_reads_120x.fasta.headerless.gz -d False
