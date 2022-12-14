import numpy as np
from sklearn.manifold import TSNE
from utils import *

@timer_func
def sort_by_tsne(infile, outfile, chunk_size):
    """Sort the reads in the input file by t-SNE.

    Args:
        infile (str): Path to the input file.
        outfile (str): Path to the output file.
        chunk_size (int): Number of reads to compute at once using t-SNE. Useful if the file can't hold in memory.
    """
    #kmer embedding sorting
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
            i, data, matrix = ts_sort(i, matrix, data, chunk_size, kmers, outfile)
    if i > 0:
        i, data, matrix = ts_sort(i, matrix, data, chunk_size, kmers, outfile)

def ts_sort(i, matrix, data, chunk_size, kmers, outfile):
    """Compute the t-SNE on the current block of reads, sort and write them to output.

    Args:
        i (int): Number of reads in the current block.
        matrix (np.array): Matrix of encoded reads.
        data (list): List of reads.
        chunk_size (int): Number of reads to compute at once using t-SNE. Useful if the file can't hold in memory.
        kmers (list): List of all kmers of length k.
        outfile (str): Path to the output file.

    Returns:
        i, data, matrix (tuple): Tuple of the number of reads in the current block, the list of reads and the matrix of encoded reads.

    """
    tsne = TSNE(n_components=1, perplexity=40, n_iter=300, init='pca')
    TS = tsne.fit_transform(matrix)
    for pc, seq, j in zip(TS, data, range(chunk_size)):
        data[j] = (pc, seq)
    data = sorted(data, key=lambda x: x[0])
    with open(outfile, 'a') as f:
        for couple in data:
            f.write(couple[1] + '\n')
    i = 0
    data = []
    matrix = np.zeros((chunk_size, len(kmers)), dtype=int)
    return i, data, matrix

if __name__ == "__main__":
    import os
    sort_by_tsne(infile='data/ecoli_100Kb_reads_10x.fasta', outfile="out_x.fasta", chunk_size=10000)
    print(monitor_gzip("out_x.fasta", 'data/headerless/ecoli_100Kb_reads_10x.fasta.headerless.gz'))
    os.remove("out_x.fasta")
    os.remove("out_x.fasta.gz")
    