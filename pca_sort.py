import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn import manifold
from utils import *
from pathlib import Path
import os


def enum_kmers(k):
    """Enumerate all possible kmers of length k.

    Args:
        k (int): [description]

    Returns:
        list (list): A lit of all possible kmers of length k.
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
    """Distribute the kmers in the list."""
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
    return [1 if kmer in sequence else 0 for kmer in kmers]


def sort_by_pca(infile: Path, outfile: Path, chunk_size):
    # kmer embedding sorting
    kmers = enum_kmers(3)
    i = 0
    data = []
    chunk_size = chunk_size if chunk_size != 0 else len(list(fasta_reader(infile)))
    matrix = np.zeros((chunk_size, len(kmers)), dtype=int)
    for read in fasta_reader(infile):
        sequence = get_sequence(read)
        matrix[i] = encode_kmers(kmers, sequence)
        i += 1
        data.append(sequence)
        if i == chunk_size:
            i, data, matrix = pc_sort(i, matrix, data, chunk_size, kmers, outfile)
    if i > 0:
        i, data, matrix = pc_sort(i, matrix, data, chunk_size, kmers, outfile)


@timer_func
def pc_sort(i, matrix, data, chunk_size, kmers, outfile):
    pca = PCA(n_components=1)
    PC = pca.fit_transform(matrix).flatten()
    for pc, seq, j in zip(PC, data, range(chunk_size)):
        data[j] = (pc, seq)
    data = sorted(data, key=lambda x: x[0])
    with open(outfile, "a") as f:
        for couple in data:
            f.write(couple[1] + "\n")
    i = 0
    data = []
    matrix = np.zeros((chunk_size, len(kmers)), dtype=int)
    return i, data, matrix


if __name__ == "__main__":
    sort_by_pca(infile='data/ecoli_100Kb_reads_40x.fasta', outfile="out_x.fasta", chunk_size=0)
    print(monitor_gzip("out_x.fasta", 'data/headerless/ecoli_100Kb_reads_40x.fasta.headerless.gz'))
    os.remove("out_x.fasta")
    os.remove("out_x.fasta.gz")
