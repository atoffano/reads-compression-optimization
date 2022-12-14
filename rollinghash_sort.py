from utils import *

@timer_func
def sort_by_minimizer(infile, outfile):
    """A function to sort a .fasta file with the rolling hash method, proposed by chatGPT.

    Args:
        infile (str): Input .fasta file
        outfile (str): Output .fasta file without headers
    """
    with open(infile, 'r') as file:
        reads = [line.strip() for line in file.readlines() if line[0] != '>'] # Remove headers if leftovers

    # define the length of the minimizer and the number of mismatches allowed
    MINIMIZER_LEN = 8 # Size of both minimizer and of rolling window
    MAX_MISMATCHES = 1

    # compute the minimizers for each read
    minimizers = []
    for read in reads:
        # initialize the rolling hash
        rolling_hash = 0
        window = read[:MINIMIZER_LEN]
        for char in window:
            rolling_hash = rolling_hash * 31 + ord(char) # ord() returns the unicode code point for a one-character string
        minimizer = rolling_hash

        # iterate over the remaining characters in the read
        for i in range(MINIMIZER_LEN, len(read)):
            # remove the leftmost character from the hash
            rolling_hash -= 31 ** (MINIMIZER_LEN - 1) * ord(read[i - MINIMIZER_LEN])
            # add the next character to the hash
            rolling_hash = rolling_hash * 31 + ord(read[i])
            # update the minimizer if necessary
            if rolling_hash < minimizer:
                minimizer = rolling_hash
        minimizers.append(minimizer)    # Sort dictionary by minimizers and generate new .fasta file with reads in new order
    sorted_reads = [x for _, x in sorted(zip(minimizers, reads))] # Sorts the reads according to the minimizers

    # write the reordered reads to a new .fasta file
    with open(outfile, 'a') as handle:
        for read in sorted_reads:
            handle.write(read + "\n")

if __name__ == '__main__':
    import os
    sort_by_minimizer('data/ecoli_100Kb_reads_40x.fasta', 'out_x.fasta')
    print(monitor_gzip("out_x.fasta", 'data/headerless/ecoli_100Kb_reads_40x.fasta.headerless.gz'))
    os.remove("out_x.fasta")
    os.remove("out_x.fasta.gz")
