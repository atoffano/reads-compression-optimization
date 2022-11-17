from typing import Tuple, Generator

from profiling import monitor
from utils import fasta_reader, get_label, get_sequence


def frequence_classification(freq_div: int = 40) -> Tuple[dict, list]:
    """Build a dict with nt frequency as key and list of reads as value, also build and
       reorder a list of those frequency

    Args:
        freq_div (int, optional): number of frequence division. Defaults to 40.

    Returns:
        Tuple[dict, list]: dictionary of the reads for each frequency and ordered list
                           of the frequencies
    """
    reader: Generator = fasta_reader("data/ecoli_100Kb_reads_120x.fasta")
    freq_dict: dict = {}
    sorted_freq: list = []

    for read in reader:
        seq: str = get_sequence(read)
        key: str = " ".join(
            [
                str(int(get_freq_interval(freq_div, seq.count(nt))))
                for nt in ["A", "T", "G", "C"]
            ]
        )

        if not key in freq_dict:
            freq_dict[key] = [read]
            sorted_freq.append(key)
        else:
            freq_dict[key].append(read)

    sorted_freq.sort()  # sort the frequence sets alphabeticaly
    return freq_dict, sorted_freq


def outfile(freq_dict: dict, sorted_freq: list) -> None:
    """write the sorted fasta file from frequency dicctionary and the sorted freq

    Args:
        freq_dict (dict): dictionary containing the reads classify by nt frequencies
        sorted_freq (list): sorted list of the frequencies found in original file
    """
    with open("data/method_1_out_120f.fasta", "w") as wipe:
        wipe.write("")  # cleaning the file

    with open("data/method_1_out_120f.fasta", "a") as file:
        for key in sorted_freq:
            for read in freq_dict[key]:
                file.write(read + "\n")


def get_freq_interval(freq_div: int, freq: float) -> float:
    """Return the

    Args:
        freq_div (int): _description_
        freq (float): _description_

    Returns:
        float: _description_
    """
    intervals: list = [(100 / freq_div) * i for i in range(freq_div)]
    prev_limit: float = 0.0
    for limit in intervals:
        if freq <= limit:
            return prev_limit
        prev_limit = limit


outfile(*frequence_classification())
