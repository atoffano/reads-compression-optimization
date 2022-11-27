# PCA sorting
import pca_sort
import method_1_jules
import method_2_antoine
import method_2_jules
from profiling import monitor
from utils import *

import glob, os
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn import manifold

if __name__ == '__main__':
    files = [
        'data/ecoli_100Kb_reads_5x.fasta',
        'data/ecoli_100Kb_reads_10x.fasta',
        'data/ecoli_100Kb_reads_20x.fasta',
        'data/ecoli_100Kb_reads_40x.fasta',
        'data/ecoli_100Kb_reads_80x.fasta',
        'data/ecoli_100Kb_reads_120x.fasta'
    ]
    log = {}
    for file in files:
        print('file= '+file)
        print('compare_to= ' + f'{file.replace("/", "/headerless/")}.headerless.gz')
        for chunk_size in range(40000, 80000, 40000):
            log_monitor_func = monitor(
                func=pca_sort.sort_by_pca(file, "out_x.fasta", chunk_size),
                input_file='out_x.fasta',
                compare_to=f'{file.replace("/", "/headerless/")}.headerless.gz'
                )
            log[(file, chunk_size)] = log_monitor_func
            os.remove('out_x.fasta')
            os.remove('out_x.fasta.gz')
    