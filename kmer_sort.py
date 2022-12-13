import random

from pathlib import Path
from typing import Generator, Tuple

from utils import *


def __get_random_kmer(size: int, infile: Path) -> str:
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


def __get_identifier(read: str, kmer: str, intervals: list) -> str:
    """return the read identifier with the occurence and position of a kmer in the read

    Args:
        read (str): the read we loof for
        kmer (str): the kmer we search in the read
        intervals (list): the position intervals used for relative position of the kmer

    Returns:
        str: the identifier computed if there is no kmer in the read it will be "0"
    """
    seq: str = get_sequence(read)
    kmer_positions: str = ""
    last_kmer_index: int = 0
    move_index: int = 1

    while move_index != 0:  # if the result of the find() is -1 we stop the loop

        # look for the next kmer starting after the previous found kmer start
        move_index = seq[last_kmer_index:].find(kmer) + 1

        if move_index > 0:
            last_kmer_index += move_index

            # concatenate str with the relative position of the kmer start
            if last_kmer_index - 1 >= intervals[-1]:
                kmer_positions += str(len(intervals))
            else:
                kmer_positions += str(
                    [j for j, i in enumerate(intervals) if last_kmer_index >= i][-1] + 1
                )

    # len("") = 0 so if there is no kmer found in the read the identifier is "0"
    return str(len(kmer_positions)) + kmer_positions


def __get_identifiers_dict(
    kmer: str,
    current_file: Path,
    next_file: Path,
    intervals_number: int,
    seqlen: int,
) -> Tuple[dict, list]:
    """Return a dictionary of the identifier and the reads associated to it, also return
       the list of unsorted reads

    Args:
        kmer (str): he kmer we search in the read
        current_file (Path): the current source file Path
        next_file (Path): the next source file Path
        intervals_number (int): number of relative position interval to use
        seqlen (int): length of the reads in the source file

    Returns:
        Tuple[dict, list]: dictionary of the identifier and there list of reads,
                           the list of unsorted reads
    """
    not_sorted: list = []
    identifiers_dict: dict = {}
    reader: Generator = fasta_reader(current_file)
    # with an interval number of 4 and seqlen of 100, intervals will look like :
    # [0.0, 25.0, 50.0, 75.0]
    intervals: list = [i * (seqlen / intervals_number) for i in range(intervals_number)]

    with open_wipe_and_add(next_file) as temp_file:
        for read in reader:
            identifier: str = __get_identifier(
                read=read, kmer=kmer, intervals=intervals
            )
            if identifier == "0":  # case of no kmer found in the read
                not_sorted.append(get_sequence(read))
                temp_file.write(read + "\n")
            else:
                if (
                    not identifier in identifiers_dict
                ):  # first occurence of the identifier
                    identifiers_dict[identifier] = [get_sequence(read)]
                else:  # already seen identifier
                    identifiers_dict[identifier].append(get_sequence(read))

    return identifiers_dict, not_sorted


def __sorting(
    infile: Path,
    output: Path,
    original_size: int,
    intervals_number: int,
    seqlen: int,
    cutoff: int,
) -> None:
    """Take the input and write the sorted output file with the given parameters

    Args:
        infile (Path): Source file Path
        output (Path): Output file Path
        original_size (int): kmer sie passed in argument
        intervals_number (int): number of relative position interval to use
        seqlen (int): length of the reads in the source file
        cutoff (int): maximum authorized kmer generation
    """
    used_kmer: list = []
    size: int = original_size
    current_file: Path = infile
    next_file: Path = Path("kmer_sort_temp_no_init")
    kmer: str = __get_random_kmer(size=size, infile=infile)

    with open_wipe_and_add(output) as output_file:
        for round in range(cutoff):

            used_kmer.append(kmer)  # keep in memory already used kmer

            # get identifier and there list of read for the current kmer, also get the
            # list of reads where no kmer were found
            identifiers_dict, not_sorted = __get_identifiers_dict(
                kmer=kmer,
                current_file=current_file,
                next_file=next_file,
                intervals_number=intervals_number,
                seqlen=seqlen,
            )

            # Write sorted reads for the current kmer in the output file
            for identifier in identifiers_dict:
                for seq in identifiers_dict[identifier]:
                    output_file.write(seq + "\n")

            # If all the reads has been sorted we end the loop
            if not not_sorted:
                break

            # Make sure we remove only temporary file and not the
            if current_file.name != infile.name:
                os.remove(current_file)

            # We put last next_file as current file and generate another next_file
            current_file = next_file
            next_file = Path(current_file.name.split("no")[0] + f"no{str(round)}")

            # Get the next kmer, if we draw an already seen kmer, we draw another with
            # a size greater of 1 to avoid (close to) infinite draw in case of small kmer
            current_size = size
            while kmer in used_kmer:
                kmer: str = __get_random_kmer(size=size, infile=current_file)
                size += 1
            size = current_size

        # if we reach the cutoff we keep the rest of the reads as an unsorted list
        if not_sorted:
            print(f"reaching {cutoff}")
            for seq in not_sorted:
                output_file.write(seq + "\n")

    # remove the lasts temporary files
    if current_file.name != infile.name:
        os.remove(current_file)
        os.remove(next_file)


@timer_func
def sort_by_kmer(
    infile: Path, output: Path, size: int, intervals_number: int, cutoff: int
) -> None:
    """Launching function of the method, check for error and sort the source file

    Args:
        infile (Path): Source file Path
        output (Path): Output file Path
        size (int): size of kmer to use
        intervals_number (int): number of relative position interval to use
        cutoff (int): maximum authorized kmer generation
    """
    random.seed(42)

    # check for error in argument values
    __error_handling(
        infile=infile,
        size=size,
        intervals_number=intervals_number,
    )

    cutoff = len(list(fasta_reader(infile))) if cutoff <= 0 else cutoff
    seqlen = len(get_sequence(next(fasta_reader(infile))))
    __sorting(
        infile=infile,
        output=output,
        original_size=size,
        intervals_number=intervals_number,
        seqlen=seqlen,
        cutoff=cutoff,
    )


def __error_handling(infile: Path, size: int, intervals_number: int) -> None:
    """Handle error than can occured from wrong arguments values

    Args:
        infile (Path): Source file Path
        size (int): kmer size passed in argument
        intervals_number (int): number of interval passed in argument
    """
    reader: Generator = fasta_reader(infile)
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
