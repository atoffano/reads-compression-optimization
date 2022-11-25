from utils import *

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

def single_kmer_sorting(file, kmer):
    outleft = []
    outfile = f'{kmer}_1.fasta'
    with open(outfile, 'a') as f:
        for read in fasta_reader(file):
            sequence = get_sequence(read)
            if kmer in sequence:
                f.write(sequence + '\n')
            else:
                outleft.append(sequence)
        for j in outleft:
            f.write(j + '\n')
    return outfile

def multiple_kmer_sorting(file):
    kmers = enum_kmers(3)
    kmer = random.choice(kmers)
    outfile = single_kmer_sorting(file, kmer)
    filelist = [outfile]
    for i in range(4):
        kmer = random.choice(kmers)
        file = f'{kmer}_{i+2}.fasta'
        outleft = []
        cnt, cnt1 = 0, 0
        with open(file, 'a') as f:
            for sequence in fasta_reader_headerless(outfile):
                if kmer in sequence:
                    cnt = cnt + 1
                    f.write(sequence + '\n')
                else:
                    cnt1 = cnt1 + 1
                    outleft.append(sequence)
            for k in outleft:
                f.write(k + '\n')
        filelist.append(file)
    print(monitor_gzip(file, 'data/headerless/ecoli_100Kb_reads_80x.fasta.headerless.gz'))
    for ent in filelist:
        print(ent)
        os.remove(ent)
    os.remove(filelist[-1] + '.gz')
    print(file)

def embed_sequence(sequence: str, kmers) -> list:
    """Embed a sequence into a list of k-mers

    Args:
        sequence (str): a sequence
        k (int): kmer size

    Returns:
        list: list of k-mers
    """
    for i in range(len(sequence)):
        embed = []
        try:
            embed.append(kmers.index(sequence[i:i+len(kmers[0])]))
        except:
            continue
    return embed

if __name__ == "__main__":
    #Single kmer sorting
    # print(monitor(single_kmer_sorting('data/ecoli_100Kb_reads_80x.fasta', 'out80x.fasta', 'ATCG')))
    # print(monitor_gzip('out80x.fasta', 'data/headerless/ecoli_100Kb_reads_80x.fasta.headerless.gz'))
    # os.remove('out80x.fasta')
    # os.remove('out80x.fasta.gz')

    #Multiple kmer sorting
    #print(monitor(multiple_kmer_sorting('data/ecoli_100Kb_reads_80x.fasta')))
    # print(monitor_gzip('out80x.fasta', 'data/headerless/ecoli_100Kb_reads_80x.fasta.headerless.gz'))

    #kmer embedding sorting
    kmers = enum_kmers(3)
    for read in fasta_reader("data/ecoli_100Kb_reads_80x.fasta"):
        sequence = get_sequence(read)
        print(list(embed_sequence(sequence, kmers)))

    