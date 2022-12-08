python read_sort.py -i data/humch1_1Mb_reads_5x.fasta  -m pca_sort -cs 50000 -c data/headerless/humch1_1Mb_reads_5x.fasta.headerless.gz --delete_output True
python read_sort.py -i data/humch1_1Mb_reads_10x.fasta  -m pca_sort -cs 100000 -c data/headerless/humch1_1Mb_reads_10x.fasta.headerless.gz --delete_output True
python read_sort.py -i data/humch1_1Mb_reads_20x.fasta  -m pca_sort -cs 200000 -c data/headerless/humch1_1Mb_reads_20x.fasta.headerless.gz --delete_output True
python read_sort.py -i data/humch1_1Mb_reads_40x.fasta  -m pca_sort -cs 400000 -c data/headerless/humch1_1Mb_reads_40x.fasta.headerless.gz --delete_output True
python read_sort.py -i data/humch1_1Mb_reads_80x.fasta  -m pca_sort -cs 800000 -c data/headerless/humch1_1Mb_reads_80x.fasta.headerless.gz --delete_output True
python read_sort.py -i data/humch1_1Mb_reads_120x.fasta  -m pca_sort -cs 1200000 -c data/headerless/humch1_1Mb_reads_120x.fasta.headerless.gz --delete_output True
python read_sort.py -i data/metagenomic_sample.fasta  -m pca_sort -cs 877179 -c data/headerless/metagenomic_sample.fasta.headerless.gz --delete_output True
