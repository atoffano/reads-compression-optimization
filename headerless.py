from utils import *

def make_headerless(infile, outfile):
    with open(infile, "r") as f:
        with open(outfile, "w") as f2:
            for line in f:
                if line.startswith(">"):
                    continue
                f2.write(line)


clean_file('data/metagenomic_sample.fasta', 'metagenomic_sample.fasta')
make_headerless('metagenomic_sample.fasta', 'data/headerless/metagenomic_sample.fasta.headerless')