from pathlib import Path


def fasta_reader(filename: Path) -> str:
    """Generator yielding read sequence with it label

    Args:
        filename (Path): path to fasta file

    Yields:
        str: labeled read
    """
    with open(filename) as file:
        lab_read: str = file.readline() + file.readline().rstrip()  # (labeled_read)
        while lab_read:
            yield lab_read
            lab_read: str = file.readline() + file.readline().rstrip()


def get_label(labeled_read: str) -> str:
    """return label of the given read

    Args:
        labeled_read (str): a labeled read

    Returns:
        str: the label of the read
    """
    return labeled_read.split("\n")[0]


def get_sequence(labeled_read: str) -> str:
    """return sequence of the given read

    Args:
        labeled_read (str): a labeled read

    Returns:
        str: the sequence of the read
    """
    return labeled_read.split("\n")[1]

############################################################################################################
# Cosine similarity sorting
############################################################################################################

def convert_to_int(str):
    match = { 'A': 0, 'T': 1, 'C': 2, 'G': 3 }
    for i in range(len(str)):
        str[i] = match[str[i]]
    return str


############################################################################################################
# Kmer sorting
############################################################################################################

def kmer_sorting(filename: Path, k: int) -> None:
    """Sort reads by kmer

    Args:
        filename (Path): path to fasta file
        k (int): kmer size
    """
    kmer_dict: dict = {}

def nb_kmers(k):
    if k < 0:
        raise ValueError()
    elif k == 0:
        return 0
    return 4**k

def enum_kmers(k):
    if k<0:
        raise ValueError()
    elif k == 0:
        return []
    list = ['']
    i = 0
    while i<k:
        list = distribution(list)
        i+=1
    return list

def distribution(list):
    adn = ['A','T','C','G']
    res = []
    i = 0
    while i < len(list):
        j = 0
        while j < 4:
            sol = list[i] + adn[j]
            res.append(sol)
            j+=1
        i+=1
    return res

