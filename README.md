# reads-compression-optimization
Optimizing read ordering for a better compression efficiency.

## Goal
New gen sequencing requires a large disk space in order to store the read output. In this project we try to reorder the reads contained in an output file to have similar reads close by.
This repository contains a few functions that utilize different strategies in order to match reads.
The functions employed here are based on two main strategies : either a kmer-based sorting strategy or a dimension reduction method.
All functions are ra

## Strategies
- pca_sort : A function that uses PCA in order to assert their similarity to each other.
- tsne_sort : Works quite like pca_sort, but with far better results and far worse execution time.
