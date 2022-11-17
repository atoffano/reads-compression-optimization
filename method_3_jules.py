from utils import fasta_reader, get_label, get_sequence
import random


def get_random_kmer():
    pass


def get_kmer_count_dict(kmer: str, file: str):
    kmer_count_dict = {}
    reader = fasta_reader(file)
    for read in reader:
        kmer_count = str(get_sequence(read).count(kmer))
        if not kmer_count in kmer_count_dict:
            kmer_count_dict[kmer_count] = [read]
        else:
            kmer_count_dict[kmer_count].append(read)

    return kmer_count_dict


def get_kmer_dict(kmer: str, file: str):

    kmer_dict = {}
    round_limit = 100
    breaking = False
    for _ in range(round_limit):
        kmer_dict[kmer] = get_kmer_count_dict(kmer, file)

        if not "0" in kmer_dict[kmer]:
            breaking = True
            break

        else:
            print(_)
            print(len(kmer_dict[kmer]["0"]))
            new_kmer_read = kmer_dict[kmer]["0"][0]
            index = random.randint(0, len(new_kmer_read) - len(kmer))
            new_kmer = new_kmer_read[index : index + len(kmer)]

            with open("method_3_temp", "w") as wipe:
                wipe.write("")
            with open("method_3_temp", "a") as file:
                for read in kmer_dict[kmer]["0"]:
                    file.write(read + "\n")

            last_zero = kmer_dict[kmer]["0"]
            del kmer_dict[kmer]["0"]

        file = "method_3_temp"
        kmer = new_kmer

    if not breaking:
        kmer_dict["LAST"] = {"0": last_zero}
    return kmer_dict


# " ": ATGC
# "b": ATGCGC
# "c": ATG

# 120 avec ATGCGC
def write_outfile():

    kmer_dict = get_kmer_dict(kmer="ATGCGC", file="data/ecoli_100Kb_reads_120x.fasta")

    with open("data/method_3_out_120.fasta", "w") as wipe:
        wipe.write("")
    with open("data/method_3_out_120.fasta", "a") as file:
        for kmer in kmer_dict:
            for count in kmer_dict[kmer]:
                for read in kmer_dict[kmer][count]:
                    file.write(read + "\n")


write_outfile()
