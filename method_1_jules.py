from typing import Tuple, Generator
from pathlib import Path
from argparse import ArgumentParser

from profiling import monitor
from utils import fasta_reader, open_wipe_and_add, get_sequence

# sort the reads by overall frequences of ACGT


def frequence_classification(input_file: Path, freq_cut: int) -> Tuple[dict, list]:
    """Build a dict with nt frequency as key and list of reads as value, also build and
       reorder a list of those frequency

    Args:
        freq_div (int, optional): number of frequence division. Defaults to 40.

    Returns:
        Tuple[dict, list]: dictionary of the reads for each frequency and ordered list
                           of the frequencies
    """
    reader: Generator = fasta_reader(filename=input_file)
    freq_dict: dict = {}
    sorted_freq: list = []

    for read in reader:
        seq: str = get_sequence(read)
        key: str = " ".join(
            [
                str(int(get_freq_interval(freq_cut, seq.count(nt))))
                for nt in ["A", "T", "G", "C"]
            ]
        )  # key is the 4 frequence interval for each nt such as '25 0 50 25'

        if not key in freq_dict:
            # if frequence set do not exist we create a list with only the read as value
            freq_dict[key] = [read]
            sorted_freq.append(key)
        else:
            # if frequence set exist we add the read to the list
            freq_dict[key].append(read)

    sorted_freq.sort()  # sort the frequence sets alphabeticaly
    return freq_dict, sorted_freq


def get_freq_interval(freq_cut: int, freq: float) -> float:
    """Return the

    Args:
        freq_cut (int): _description_
        freq (float): _description_

    Returns:
        float: _description_
    """
    intervals: list = [(100 / freq_cut) * i for i in range(freq_cut)]
    prev_limit: float = 0.0
    for limit in intervals:
        if freq <= limit:
            return prev_limit
        prev_limit = limit


def write_outfile(in_file: Path, out_file: Path, freq_cut: Path) -> None:
    """write the sorted fasta file from frequency dicctionary and the sorted freq

    Args:
        freq_dict (dict): dictionary containing the reads classify by nt frequencies
        sorted_freq (list): sorted list of the frequencies found in original file
    """
    freq_dict, sorted_freq_list = frequence_classification(
        input_file=in_file, freq_cut=freq_cut
    )

    with open_wipe_and_add(filename=out_file) as file:
        # group are write in the file folowing the ordered list
        for key in sorted_freq_list:
            for read in freq_dict[key]:
                file.write(read + "\n")


def main():
    parser = ArgumentParser(
        prog="method_2_jules.py",
        description="",
        epilog="Text at the bottom of help",
    )

    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    parser.add_argument("-c", "--freq_cut")

    args = parser.parse_args()

    freq_cut: int = 40 if args.freq_cut == None else int(args.freq_cut)

    print(f"Parameters are :\nFrequence cutting : {freq_cut}")

    write_outfile(in_file=args.input, out_file=args.output, freq_cut=freq_cut)


if __name__ == "__main__":
    main()
