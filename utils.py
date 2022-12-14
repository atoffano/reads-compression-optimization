import time
import os, gzip
import shutil

from pathlib import Path
from io import TextIOWrapper

def timer_func(func):
    """Decorator meant to be used on functions to show the execution time of the function.
    This is our primary way of measuring execution time.

    Args:
        func (Function): The function to be timed

    Returns:
        wrap_func: The function wrapped insinde our decorator.
    """
    # This function shows the execution time of
    # the function object passed
    def wrap_func(*args, **kwargs):
        t1 = time.time()
        result = func(*args, **kwargs)
        t2 = time.time()
        print(f"Function {func.__name__!r} executed in {(t2-t1):.4f}s")
        return result

    return wrap_func


def monitor_gzip(input_file, compare_to):
    """Function to monitor the compression ratio of the output file.
    
    Args:
        input_file (Path): path to the output file
        compare_to (Path): path to the input file
    """
    with open(input_file, "rb") as f_in:
        with gzip.open(f"{input_file}.gz", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    ret = os.path.getsize(compare_to) / os.path.getsize(f"{input_file}.gz")
    return ret


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
    """Function to remove headers from ."""
    with open(filein, "r") as f:
        with open(fileout, "w") as f2:
            for line in f:
                if line.startswith(">"):
                    continue
                f2.write(line)

def clean_file(filein, fileout):
    """Function to clean the 'mix' testing file and output it in a correct format."""
    with open(filein) as fi:
        with open_wipe_and_add(fileout) as fo:
            line = fi.readline()
            seq = ""
            while line:
                if line.startswith(">"):
                    if seq != "":
                        fo.write("".join(seq.split("\n")) + "\n")
                    seq = ""
                    fo.write(line)
                else:
                    seq += line
                line = fi.readline()
            fo.write(seq)


def open_wipe_and_add(filename: Path) -> TextIOWrapper:
    """Clean the file and open it as "a"
    Args:
        filename (Path): the output_file
    Returns:
        TextIOWrapper: the opened file
    """
    with open(filename, "w") as wipe:
        wipe.write("")

    return open(filename, "a")

def clean_file(filein, fileout):
    with open(filein) as fi:
        with open_wipe_and_add(fileout) as fo:
            line = fi.readline()
            seq = ""
            while line:
                if line.startswith(">"):
                    if seq != "":
                        fo.write("".join(seq.split("\n")) + "\n")
                    seq = ""
                    fo.write(line)
                else:
                    seq += line
                line = fi.readline()
            fo.write(seq)

# Functions for dimension reduction methods.  
def enum_kmers(k):
    """Enumerate all possible kmers of length k.

    Args:
        k (int): [description]

    Returns:
        list (list): A list of all possible kmers of length k.
    """
    if k<0:
        raise ValueError()
    elif k == 0:
        return []
    list = [""]
    i = 0
    while i < k:
        list = distribution(list)
        i += 1
    return list

def distribution(list):
    """Distribute the kmers in the list.
    
    Args:
        list (list): A list of kmers.

    Returns:
        list (list): A list of kmers with all possible kmers of length k.
    """
    adn = ['A','T','C','G']
    res = []
    i = 0
    while i < len(list):
        j = 0
        while j < 4:
            sol = list[i] + adn[j]
            res.append(sol)
            j += 1
        i += 1
    return res

def encode_kmers(kmers, sequence):
    """Encode the kmers in the sequence.
        
    Args:
        kmers (list): List of all kmers of length k.
        sequence (str): Sequence to encode.

    Returns:
        list (list): A list of 1s and 0s, 1 if the kmer is in the sequence, 0 otherwise.
    """
    return [1 if kmer in sequence else 0 for kmer in kmers]


if __name__ == "__main__":
    pass