
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn import manifold
from utils import *
import os

def enum_kmers(k):
    if k<0:
        raise ValueError()
    elif k == 0:
        return []
    list = ['']
    i = 0
    while i<k:
        list = distribution(list)
        i+=1
    return list

def distribution(list):
    adn = ['A','T','C','G']
    res = []
    i = 0
    while i < len(list):
        j = 0
        while j < 4:
            sol = list[i] + adn[j]
            res.append(sol)
            j+=1
        i+=1
    return res

def encode_kmers(kmers, sequence):
    return [1 if kmer in sequence else 0 for kmer in kmers]

def encode_to_kmers(kmers, sequence):
    encoding = []
    for i in range(len(sequence) - len(kmers[0]) + 1):
        kmer = sequence[i:i+len(kmers[0])]
        if kmer in kmers:
            encoding.append(kmers.index(kmer))
    return encoding

def sort_by_pca(infile, outfile, chunk_size):
    #kmer embedding sorting
    kmers = enum_kmers(3)
    i = 0
    data = []
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
    if not os.path.exists(outfile):
        print('chunk_size of size ' + str(chunk_size) + 'is larger than the number of reads. Using chunk_size equal to the nb of reads ('+str(i)+') instead')
        i, data, matrix = pc_sort(i, matrix, data, chunk_size, kmers, outfile)


def pc_sort(i, matrix, data, chunk_size, kmers, outfile):
    pca = PCA(n_components=1)
    PC = pca.fit_transform(matrix).flatten()
    for pc, seq, j in zip(PC, data, range(chunk_size)):
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
    sort_by_pca(infile='data/ecoli_100Kb_reads_120x.fasta', outfile="out_x.fasta", chunk_size=12000)
    print(monitor_gzip("out_x.fasta", 'data/headerless/ecoli_100Kb_reads_120x.fasta.headerless.gz'))

    
