from pathlib import Path
from profiling import monitor, monitor_gzip
import os
import random
import glob

def fasta_reader(filename: Path) -> str:
    """Generator yielding read sequence with it label

    Args:
        filename (Path): path to fasta file

    Yields:
        str: labeled read
    """
    with open(filename) as file:
        lab_read: str = file.readline() + file.readline().rstrip()  # (labeled_read)
        while lab_read:
            yield lab_read
            lab_read: str = file.readline() + file.readline().rstrip()


def get_label(labeled_read: str) -> str:
    """return label of the given read

    Args:
        labeled_read (str): a labeled read

    Returns:
        str: the label of the read
    """
    return labeled_read.split("\n")[0]


def get_sequence(labeled_read: str) -> str:
    """return sequence of the given read

    Args:
        labeled_read (str): a labeled read

    Returns:
        str: the sequence of the read
    """
    return labeled_read.split("\n")[1]

def fasta_reader_headerless(filename: Path) -> str:
    """Generator yielding read sequence from headerless file

    Args:
        filename (Path): path to fasta file

    Yields:
        str: labeled read
    """
    with open(filename) as file:
        read: str = file.readline().rstrip()  # (labeled_read)
        while read:
            yield read
            read: str = file.readline().rstrip()

def remove_headers(filein, fileout):
    with open(filein, 'r') as f:
        with open(fileout, 'w') as f2:
            for line in f:
                if line.startswith('>'):
                    continue
                f2.write(line)

############################################################################################################
# Cosine similarity sorting
############################################################################################################

def convert_to_int(str):
    match = { 'A': 0, 'T': 1, 'C': 2, 'G': 3 }
    for i in range(len(str)):
        str[i] = match[str[i]]
    return str


