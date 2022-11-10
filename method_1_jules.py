from utils import fasta_reader, get_label, get_sequence
from pathlib import Path


# 0 ATAATTATA <- bot
# 1 ATTAAATAT
# 2 TATACTCCT
# 3 TATCGCATC <- mid (ou mid_bot / mid_top si impair)
# 4 GCGGACCAG
# 5 GGAGGAGGG
# 6 GGAGGGGAG <- top

# ici imaginons que mid > bot > top :
#
#       on place entre 2 et 3
#
#               OU
#
#       on refait avec seulement les seqeunce 0,1,2,3 pour plus préçis
#
# on peut spécifier en paramètres un nombre de "découpe" par séquence
#
# 2e découpe ici si mid_bot > mid_top > bot > top
# alors on le met entre 1 et 2

# on fait 3 ou 4 * nombre de découpe comparaison par séquence


def get_comparisons():
    reader = fasta_reader("ecoli_100Kb_reads_5x.fasta")


def get_indexes(reader):
    mid = None
    mid_bot = None
    mid_top = None
    i = 0
    for read in reader:
        i += 


def compare_sequence():
    pass


def insert_sequence():
    pass
