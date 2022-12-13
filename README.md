# reads-compression-optimization
Optimizing read ordering for a better compression efficiency.

## Goal
Sequencing data requires a large disk space in order to store its output. In this project we try to reorder the reads contained in a sequencing output file in a .fasta format to have similar reads close by.

This repository contains a few functions that utilize different strategies in order to match reads.
The functions employed here are based on two main strategies : either a minimizer-based sorting method or a dimension reduction method.

## Strategies
- `kmer_sort` : A method that compresses the file by matching reads through a k-mer minimizer.
- `rollinghash_sort` : Uses a rolling hash function in order to find the best minimizer to reorder the reads.
- `pca_sort` : A function that uses PCA in order to assert their similarity to each other.
- `tsne_sort` : Works quite like pca_sort, but with far better results and far worse execution time.

|             |    Kmer    | Pollinghash |    PCA    |   t-SNE   |
|-------------|------------|-------------|-----------|-----------|
| Compression | ●●         | ●●●●        |   ●●      | ●●●●      |
|    Speed    | ●●         | ●●●●        |   ●       | ●●        |
|  RAM Usage  | ●●         | ●●●●        |   ●●●     | ●         |

## Usage
`$ python read_sort.py --input [input fasta] --method [strategy_name]`

### Arguments
 --compare_to [.fasta.gz of the input without headers] 
 --chunk_size : 
 --delete_output [True/False] : Whether to delete the output files.
