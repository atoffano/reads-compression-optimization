import numpy as np
import random
from utils import *
import os


def hammingDist(str1, str2):
    i = 0
    count = 0
    while(i < len(str1)):
        if(str1[i] != str2[i]):
            count += 1
        i += 1
    return count

def update_ref(match, pool_ratios, pool_size):
    ref = ''
    for i in range(len(match)):
        pool_ratios[i][match[i]] += 1 / pool_size
        for key, val in pool_ratios[i].items():
            if val == max(pool_ratios[i].values()):
                ref += key
                continue
    pool_size += 1
    return ref, pool_ratios, pool_size

def complementary_strand(strand):
    comp = ''
    for i in strand:
        if i == 'A':
            comp += 'T'
        elif i == 'T':
            comp += 'A'
        elif i == 'C':
            comp += 'G'
        elif i == 'G':
            comp += 'C'
    return comp

# @timer_func
# def get_ref_read(infile):
#     pool_ratios = {i: {'A': 0, 'T': 0, 'C': 0, 'G': 0} for i in range(100)}
#     pool_size = 1
#     ref = ''
#     for read in fasta_reader(infile):
#         sequence = get_sequence(read)
#         ref, pool_ratios, pool_size = update_ref(sequence, pool_ratios, pool_size)
#     return ref, pool_ratios, pool_size


@timer_func
def sort_by_hash(infile, outfile):
    h1= {}
    h0 = []
    for read in fasta_reader(infile):
        sequence = get_sequence(read)
        h0.append(sequence)
        # if sequence in h1:
        #     h1[sequence].append(sequence)
        # else:
        #     h1[sequence] = [sequence]
    with open(outfile, 'w') as f:
        for r in sorted(h0):
            f.write(r + '\n')

    # with open(outfile, 'w') as f:
    #     keys = h1.keys()
    #     sorted_keys = sorted(keys)
    #     for key in sorted_keys:
    #         for seq in h1[key]:
    #             f.write(seq + '\n')

# @timer_func
# def sort_by_hash(infile, outfile, maxshift=50):
#     h1, h2 = {}, {}
#     h0 = []
#     for read in fasta_reader(infile):
#         sequence = get_sequence(read)

#         h0.append(sequence)

#         if sequence[18:50] in h1:
#             h1[sequence[18:50]].append(sequence)
#         else:
#             h1[sequence[18:50]] = [sequence]

#         if sequence[50:82] in h2:
#             h2[sequence[50:82]].append(sequence)
#         else:
#             h2[sequence[50:82]] = [sequence]

#     ref = random.choice(h0) # Pick a random read as reference.
#     pool_ratios = {i: {'A': 0, 'T': 0, 'C': 0, 'G': 0} for i in range(100)}
#     pool_size = 1

#     with open(outfile, 'w') as f:
#         while len(h0) > 0:
#             for i in range(maxshift):
#                 try:
#                     S = set(h1[ref[18+i:50+i]])
#                     S.update(h2[ref[50+i:82+i]])
#                     S = list(S)
#                 except KeyError:
#                     S = []
#                 if len(S) == 0:
#                     ref = random.choice(h0)
#                     break
#                 D = random.choice(S)
#                 if hammingDist(ref[i:100], D[0:100-i]) <= 2:
#                     if D in h1[ref[18+i:50+i]]:
#                         h1[ref[18+i:50+i]].remove(D)
#                     if D in h2[ref[50+i:82+i]]:
#                         h2[ref[50+i:82+i]].remove(D)
#                     h0.remove(D)
#                     ref, pool_ratios, pool_size = update_ref(D, pool_ratios, pool_size)
#                     f.write(D + '\n')
#                     break
#             if i == maxshift:
#                 pool_ratios = {i: {'A': 0, 'T': 0, 'C': 0, 'G': 0} for i in range(100)}
#                 pool_size = 1
#                 for read in h1[ref[18:50]]:
#                     if read in h0:
#                         h0.remove(read) 
#                         f.write(read + '\n')
#                 for read in h2[ref[50:82]]:
#                     if read in h0:
#                         h0.remove(read)
#                         f.write(read + '\n')
#                 for i in range(maxshift):
#                     del h1[ref[18+i:50+i]]
#                     del h2[ref[50+i:82+i]]
#                 try:
#                     ref = random.choice(h0)
#                 except:
#                     continue


# @timer_func
# def sort_by_hash(infile, maxshift=50):
#     h1, h2 = {}, {}
#     h0 = []
#     for read in fasta_reader(infile):
#         sequence = get_sequence(read)

#         h0.append(sequence)

#         if sequence[35:50] in h1:
#             h1[sequence[45:50]].append(sequence)
#         else:
#             h1[sequence[45:50]] = [sequence]

#         if sequence[50:55] in h2:
#             h2[sequence[50:55]].append(sequence)
#         else:
#             h2[sequence[50:55]] = [sequence]



#     ref = random.choice(h0) # Pick a random read as reference.
#     pool_ratios = {i: {'A': 0, 'T': 0, 'C': 0, 'G': 0} for i in range(100)}
#     pool_size = 1
#     refr = complementary_strand(ref)
#     rev = []
#     tot_nb_reads = len(h0)

#     while len(h0) > 0:
#         for i in range(maxshift):
#             S = set(h1[ref[45+i:50+i]])
#             S.update(h2[ref[50+i:55+i]])
        
#             D = random.choice(list(S))
#             if hammingDist(ref[i:100], D[0:100-i]) <= 2:
#                 ref, pool_ratios, pool_size = update_ref(D, pool_ratios, pool_size)
#                 refr = complementary_strand(ref)
#                 rev.append(0)
#                 if D in h0:
#                     h0.remove(D)
#                 break
#             else:
#                 SR = set(h1[refr[45-i:50-i]])
#                 SR.update(h2[refr[50-i:55-i]])

#             if hammingDist(refr[0:100-i], D[i:100]) <= 2:
#                 ref, pool_ratios, pool_size = update_ref(complementary_strand(D), pool_ratios, pool_size)
#                 refr = complementary_strand(ref)
#                 rev.append(1)

#         if tot_nb_reads == len(h0):
#             ref = random.choice(list(h0.values()))
#             rev.append(0)
if __name__ == '__main__':
    sort_by_hash('data/ecoli_100Kb_reads_40x.fasta', 'out_x.fasta')
    print(monitor_gzip("out_x.fasta", 'data/headerless/ecoli_100Kb_reads_40x.fasta.headerless.gz'))
    os.remove("out_x.fasta")
    os.remove("out_x.fasta.gz")