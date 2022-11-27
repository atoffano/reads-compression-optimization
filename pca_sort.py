
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn import manifold
from utils import *
from profiling import monitor

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

def sort_by_pca(infile, outfile, chunk_size):
    #kmer embedding sorting
    kmers = enum_kmers(3)
    i = 0
    data = []
    matrix = np.zeros((chunk_size, len(kmers)), dtype=int)
    for read_nb, read in enumerate(fasta_reader(infile)):
        sequence = get_sequence(read)
        matrix[i] = [1 if kmer in sequence else 0 for kmer in kmers]
        i += 1
        data.append(sequence)
        if i == chunk_size:
            i, data, matrix = pca_sort(i, matrix, data, chunk_size, kmers, outfile)


def pca_sort(i, matrix, data, chunk_size, kmers, outfile):
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
    pass
    # print('f')
    # # print(monitor(sort_by_pca("data/ecoli_100Kb_reads_80x.fasta", "out_80x.fasta", 45000)))
    # # print(monitor_gzip('out_80x.fasta', 'data/headerless/ecoli_100Kb_reads_80x.fasta.headerless.gz'))
    # # os.remove('out.fasta')
    # # os.remove('out.fasta.gz')
    # files = [
    #     'data/ecoli_100Kb_reads_5x.fasta',
    #     'data/ecoli_100Kb_reads_10x.fasta',
    #     'data/ecoli_100Kb_reads_20x.fasta',
    #     'data/ecoli_100Kb_reads_40x.fasta',
    #     'data/ecoli_100Kb_reads_80x.fasta',
    #     'data/ecoli_100Kb_reads_120x.fasta'
    # ]

    # log = {}
    # for file in files:
    #     for chunk_size in range(40000, 80000, 40000):
    #         log_monitor_func = monitor(sort_by_pca(file, "out_x.fasta", chunk_size))
    #         log_monitor_gzip = monitor_gzip('out_x.fasta', f'{file.replace("/", "/headerless/")}.headerless.gz')
    #         log[(file, chunk_size)] = [log_monitor_func, log_monitor_gzip]
    #         print(log)
    #         os.remove('out_x.fasta')
    #         os.remove('out_x.fasta.gz')
    #         break
 
