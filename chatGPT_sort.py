import numpy as np
import random
from utils import *
import os

@timer_func
def sort_by_minimizer(infile, outfile):
    with open(infile, 'r') as file:
        reads = [line.strip() for line in file.readlines() if line[0] != '>']

    # define the length of the minimizer and the number of mismatches allowed
    MINIMIZER_LEN = 8
    MAX_MISMATCHES = 1

    # compute the minimizers for each read
    minimizers = []
    for read in reads:
        # initialize the rolling hash
        rolling_hash = 0
        window = read[:MINIMIZER_LEN]
        for char in window:
            rolling_hash = rolling_hash * 31 + ord(char)
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
    sorted_reads = [x for _, x in sorted(zip(minimizers, reads))]

    # write the reordered reads to a new .fasta file
    with open(outfile, 'a') as handle:
        for read in sorted_reads:
            handle.write(read + "\n")


@timer_func
def sort_by_minimizer2(infile, outfile):
    with open(infile, 'r') as file:
        reads = [line.strip() for line in file.readlines() if line[0] != '>']
    
        # Create a dictionary to store the minimizers and their corresponding sequences
        min_dict = {}
        
        # Iterate through the reads and find the minimizer for each sequence
        for read in reads:
            # Get the minimizer for the current sequence
            min_seq = minimizer(read)
            
            # Add the minimizer and corresponding sequence to the dictionary
            min_dict[min_seq] = read
        
        # Open a new .fasta file to write the reordered sequences
        with open(outfile, 'a') as f_out:
            # Sort the dictionary keys (minimizers) in ascending order
            sorted_keys = sorted(min_dict.keys())
            
            # Iterate through the sorted keys and write the corresponding sequences to the new .fasta file
            for key in sorted_keys:
                f_out.write(min_dict[key] + '\n')

def minimizer(sequence):
    # Set the k-mer length and window size
    k = 4
    w = 10

    # Initialize the minimizer to be the first k-mer in the sequence
    min_seq = sequence[:k]

    # Iterate through the sequence with the specified window size
    for i in range(len(sequence) - w):
        # Get the current window of the sequence
        window = sequence[i:i+w]
        
        # Iterate through the window and find the minimum k-mer
        for j in range(len(window) - k):
            curr_seq = window[j:j+k]
            if curr_seq < min_seq:
                min_seq = curr_seq

    # Return the minimizer sequence
    return min_seq


if __name__ == '__main__':
    sort_by_minimizer('data/ecoli_100Kb_reads_40x.fasta', 'out_x.fasta')
    print(monitor_gzip("out_x.fasta", 'data/headerless/ecoli_100Kb_reads_40x.fasta.headerless.gz'))
    # os.remove("out_x.fasta")
    # os.remove("out_x.fasta.gz")
