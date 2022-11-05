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
