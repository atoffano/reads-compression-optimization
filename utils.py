import time
import psutil
import multiprocessing as mp
import os, gzip
import shutil

from pathlib import Path
from io import TextIOWrapper


def timer_func(func):
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
    unsorted_comp_size = {
        "data/headerless/ecoli_100Kb_reads_10x.fasta.headerless.gz": 292244,
        "data/headerless/ecoli_100Kb_reads_120x.fasta.headerless.gz": 3508076,
        "data/headerless/ecoli_100Kb_reads_20x.fasta.headerless.gz": 586108,
        "data/headerless/ecoli_100Kb_reads_40x.fasta.headerless.gz": 1170066,
        "data/headerless/ecoli_100Kb_reads_5x.fasta.headerless.gz": 146953,
        "data/headerless/ecoli_100Kb_reads_80x.fasta.headerless.gz": 2342520,
    }
    with open(input_file, "rb") as f_in:
        with gzip.open(f"{input_file}.gz", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    ret = unsorted_comp_size[compare_to] / os.path.getsize(f"{input_file}.gz")
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
    with open(filein, "r") as f:
        with open(fileout, "w") as f2:
            for line in f:
                if line.startswith(">"):
                    continue
                f2.write(line)


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


if __name__ == "__main__":
    from tsne_sort import sort_by_tsne

    sort_by_tsne(
        infile="data/ecoli_100Kb_reads_80x.fasta",
        outfile="out_x.fasta",
        chunk_size=80000,
    )
    print(
        monitor_gzip(
            "out_x.fasta", "data/headerless/ecoli_100Kb_reads_80x.fasta.headerless.gz"
        )
    )
