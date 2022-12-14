# reads-compression-optimization
Optimizing read ordering for a better compression efficiency.

## Purpose
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
| Compression | ●●●        | ●●●●        |   ●●      | ●●●●      |
|    Speed    | ●●         | ●●●         |   ●●●●    | ●         |
|  RAM Usage  | ●●●●       | ●●●●        |   ●●      | ●         |

Each method is ranked between ● (worst) and ●●●● (best) in different metrics.

## Basic Usage
`$ python read_organizer.py --input [input fasta] --method [strategy_name]`

### Arguments
`read_organizer.py [-h] -i INPUT [-o OUTPUT] [-m METHOD] [-d DELETE_OUTPUT] [-c COMPARE_TO] [-s SIZE_KMER] [-cs CHUNK_SIZE] [-cu CUTOFF] [-in INTERVALS_NUMBER]`
```
 -h, --help            
                        show this help message and exit
  -i INPUT, --input INPUT
                        Input file path
  -o OUTPUT, --output OUTPUT
                        Output file path, if not provided output will be save in the script file folder
  -d DELETE_OUTPUT, --delete_output DELETE_OUTPUT
                        passed 'T' (True) to delete any output after the algorithm. This only should be used to monitor performance
  -m METHOD, --method METHOD
                        Method to use between 'tsne_sort', 'pca_sort','kmer_sort','rollinghash_sort'
  -c COMPARE_TO, --compare_to COMPARE_TO
                        File path to compare the compressed output file with
  -s SIZE_KMER, --size_kmer SIZE_KMER
                        Kmer size to use in kmer_sort method
  -cs CHUNK_SIZE, --chunk_size CHUNK_SIZE
                        Chunk size to use in pca_sort and tsne_sort method, by default chunk sizewill be equal to number of read in the file
  -cu CUTOFF, --cutoff CUTOFF
                        Cutoff to use in kmer_sort method, by default there is no cutoff
  -in INTERVALS_NUMBER, --intervals_number INTERVALS_NUMBER
                        Number of intervals to use in the kmer_sort method
```
