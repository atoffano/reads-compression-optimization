from pathlib import Path
from io import TextIOWrapper
import subprocess
from profiling import monitor, monitor_gzip
import os
import random
import glob


def format_output(process):
    out, error = process.communicate()
    output = str(out)
    print(output if output != "b''" else "")
    print(error if error != None else "")


def gzip_in_and_out(in_file, out_file):
    copied_in_file = "_copy_".join(in_file.split("."))
    bashCommand = f"/bin/cp {in_file} {copied_in_file}"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    format_output(process)

    for file in [copied_in_file, out_file]:
        bashCommand = f"/bin/gzip -f {file}"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        print("compressing...")
        format_output(process)

    return copied_in_file


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


############################################################################################################
# Cosine similarity sorting
############################################################################################################


def convert_to_int(str):
    match = {"A": 0, "T": 1, "C": 2, "G": 3}
    for i in range(len(str)):
        str[i] = match[str[i]]
    return str


############################################################################################################
# Kmer sorting
############################################################################################################


def kmer_sorting(filename: Path, k: int) -> None:
    """Sort reads by kmer

    Args:
        filename (Path): path to fasta file
        k (int): kmer size
    """
    kmer_dict: dict = {}


def nb_kmers(k):
    if k < 0:
        raise ValueError()
    elif k == 0:
        return 0
    return 4**k


def enum_kmers(k):
    if k < 0:
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
    adn = ["A", "T", "C", "G"]
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


def single_kmer_sorting(file, kmer):
    outleft = []
    outfile = f"{kmer}_1.fasta"
    with open(outfile, "a") as f:
        for read in fasta_reader(file):
            sequence = get_sequence(read)
            if kmer in sequence:
                f.write(sequence + "\n")
            else:
                outleft.append(sequence)
        for j in outleft:
            f.write(j + "\n")
    return outfile


def multiple_kmer_sorting(file):
    kmers = enum_kmers(3)
    kmer = random.choice(kmers)
    outfile = single_kmer_sorting(file, kmer)
    filelist = [outfile]
    for i in range(4):
        kmer = random.choice(kmers)
        file = f"{kmer}_{i+2}.fasta"
        outleft = []
        cnt, cnt1 = 0, 0
        with open(file, "a") as f:
            for sequence in fasta_reader_headerless(outfile):
                if kmer in sequence:
                    cnt = cnt + 1
                    f.write(sequence + "\n")
                else:
                    cnt1 = cnt1 + 1
                    outleft.append(sequence)
            for k in outleft:
                f.write(k + "\n")
        filelist.append(file)
    print(
        monitor_gzip(file, "data/headerless/ecoli_100Kb_reads_80x.fasta.headerless.gz")
    )
    for ent in filelist:
        print(ent)
        os.remove(ent)
    os.remove(filelist[-1] + ".gz")
    print(file)


if __name__ == "__main__":
    # print(monitor(single_kmer_sorting('data/ecoli_100Kb_reads_80x.fasta', 'out80x.fasta', 'ATCG')))
    # print(monitor_gzip('out80x.fasta', 'data/headerless/ecoli_100Kb_reads_80x.fasta.headerless.gz'))
    # os.remove('out80x.fasta')
    # os.remove('out80x.fasta.gz')

    print(monitor(multiple_kmer_sorting("data/ecoli_100Kb_reads_80x.fasta")))
    # print(monitor_gzip('out80x.fasta', 'data/headerless/ecoli_100Kb_reads_80x.fasta.headerless.gz'))
