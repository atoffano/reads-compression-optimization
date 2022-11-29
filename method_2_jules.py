import random
import subprocess
from pathlib import Path
from argparse import ArgumentParser
from typing import Generator
import subprocess
from utils import fasta_reader, get_sequence, gzip_out, open_wipe_and_add
from profiling import monitor, monitor_gzip

# sort the reads by kmers findings in it

# trouvé une méthode pour que la sélection des kmers soit smart et pas aléatoire


def get_first_kmer(kmer: str, size: str, in_file: Path) -> str:

    if kmer != None:
        return kmer  # if a kmer is give as parameter use it as first kmer
    else:
        try:
            size = int(size)
        except ValueError:
            raise ValueError("kmer size (-s, --size-mer) must be an integer")

        # if no kmer given in args, select a random kmer in the first read
        reader: Generator = fasta_reader(in_file)
        read_seq: str = get_sequence(next(reader))
        kmer_start: int = random.randint(0, len(read_seq) - size)

        return read_seq[kmer_start : kmer_start + size]


def get_kmer_count_dict(kmer: str, file: Path) -> dict:
    """Build the dict with the num of occurence of the given kmer as key and the lists
       of the correspondings reads as value

    Args:
        kmer (str): the given kmer
        file (Path): the input file path

    Returns:
        dict: dict with num of occurence as keys and reads lists as values
    """
    kmer_count_dict: dict = {}
    reader: Generator = fasta_reader(file)
    for read in reader:
        kmer_count: str = str(get_sequence(read).count(kmer))
        if not kmer_count in kmer_count_dict:
            kmer_count_dict[kmer_count] = [read]
        else:
            kmer_count_dict[kmer_count].append(read)

    return kmer_count_dict


def get_kmer_dict(kmer: str, in_file: Path) -> dict:

    kmer_dict: dict = {}
    round_limit: int = 100
    breaking: bool = False

    # a dict is build for the first kmer and then we pick another kmer to sort the
    # reads in which the previous kmer had 0 occurences
    for _ in range(round_limit):

        # build the occurence dict for the current kmer
        kmer_dict[kmer]: dict = get_kmer_count_dict(kmer, in_file)

        # if all the reads has been sorted we end the loop
        if not "0" in kmer_dict[kmer]:
            breaking: bool = True
            break

        else:
            # pick a random kmer in the first sequence of the reads that had 0
            # occurences of the previous kmer
            new_kmer_read: str = get_sequence(kmer_dict[kmer]["0"][0])
            index: int = random.randint(0, len(new_kmer_read) - len(kmer))
            new_kmer: str = new_kmer_read[index : index + len(kmer)]

            # write the reads with 0 occurence of the previous kmer as new input file
            with open("method_2_temp", "w") as wipe:
                wipe.write("")
            with open("method_2_temp", "a") as file:
                for read in kmer_dict[kmer]["0"]:
                    file.write(read + "\n")

            last_zero: list = kmer_dict[kmer]["0"]
            del kmer_dict[kmer]["0"]

        in_file: Path = "method_2_temp"
        kmer: str = new_kmer

    if not breaking:
        # if we reach 100 kmer we keep the rest of the reads as an unsorted list
        kmer_dict["LAST"]: dict = {"0": last_zero}
    return kmer_dict


def write_outfile(in_file: str, out_file: str, kmer: str = "", size: str = "") -> None:
    in_file = Path(in_file)
    out_file = Path(out_file)

    first_kmer: str = get_first_kmer(kmer=kmer, size=size, in_file=in_file)
    kmer_dict: dict = get_kmer_dict(first_kmer, in_file)
    kmer_list = [kmer for kmer in kmer_dict.keys()]
    kmer_list.sort()
    sorted_kmer_dict = {}

    for kmer in kmer_list:
        sorted_kmer_dict[kmer] = kmer_dict[kmer]

    with open_wipe_and_add(out_file) as file:
        for kmer in sorted_kmer_dict:
            for count in sorted_kmer_dict[kmer]:
                for read in sorted_kmer_dict[kmer][count]:
                    file.write(get_sequence(read) + "\n")


def erro_handling(input_file: Path, kmer: str, size: int) -> None:
    """Handle error regarding the kmer or size given in arguments

    Args:
        input_file (Path): input file
        kmer (str): given kmer
        size (int): given size

    Raises:
        ValueError: kmer size check
        ValueError: positive and non-null size check
        ValueError: size limit check
    """
    reader: Generator = fasta_reader(filename=input_file)
    read = next(reader)

    if kmer != "random" and len(kmer) >= len(get_sequence(read)):
        raise ValueError("Kmer cannot be longer than reads size")

    if size <= 0:
        raise ValueError("Size cannot be equal or lesser than 0")

    if size > len(get_sequence(read)):
        raise ValueError("Size cannot be greater than reads size")


def launch_method_2(
    input: Path,
    output: Path,
    delete_output: bool = False,
    kmer: str = "random",
    size: int = 6,
) -> dict:
    """function for launching the method inside a python script

    Args:
        input (Path): _description_
        output (Path): _description_
        delete_output (bool, optional): _description_. Defaults to False.
        kmer (str, optional): _description_. Defaults to "random".
        size (int, optional): _description_. Defaults to 6.

    Returns:
        dict: logs of the operation
    """
    erro_handling(input_file=input, kmer=kmer, size=size)

    # Reorder and get monitoring values
    monitoring_values = monitor(
        write_outfile(
            in_file=input,
            out_file=output,
            kmer=kmer,
            size=size,
        )
    )

    # Compress the files
    base_file = "data/ori/" + input.split("data/")[1] + ".gz"
    out_gz = gzip_out(output)

    # get compression rate
    monitoring_values["rate"] = monitor_gzip(out_gz, base_file)

    # Delete the file if -d is on
    if delete_output:
        bashCommand = f"rm -rf {out_gz}"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        out, error = process.communicate()

    monitoring_values["file"] = input
    monitoring_values["kmer-size"] = size
    return monitoring_values


def main():
    parser = ArgumentParser(
        prog="method_2_jules.py",
        description="",
        epilog="Text at the bottom of help",
    )

    parser.add_argument("-i", "--input")
    parser.add_argument("-d", "--delete_output", default=False)
    parser.add_argument("-o", "--output")
    parser.add_argument("-k", "--kmer", default="random")
    parser.add_argument("-s", "--size_kmer", default=6)

    args = parser.parse_args()
    size = int(args.size_kmer)
    delete_output = False if args.delete_output == "False" else True
    # Check for errors
    erro_handling(input_file=args.input, kmer=args.kmer, size=size)

    # Print parameters
    print(f"Parameters are :\n\tFirst kmer : {args.kmer}\n\tKmer size : {size}\n")

    # Reorder and get monitoring values
    print("reordering...\n")
    monitoring_values = monitor(
        write_outfile(
            in_file=args.input,
            out_file=args.output,
            kmer=args.kmer,
            size=size,
        )
    )

    # Print monitoring values
    for key in monitoring_values:
        print(f"{key} : {monitoring_values[key]}")

    # Compress the files
    base_file = "data/ori/" + args.input.split("data/")[1] + ".gz"
    out_gz = gzip_out(args.output)

    # Print compression rate
    compression_rate = monitor_gzip(out_gz, base_file)

    # Delete the file if -d is on
    if args.delete_output:
        bashCommand = f"rm -rf {out_gz}"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        out, error = process.communicate()

    with open_wipe_and_add("log_method2.txt") as file:
        for key in monitoring_values:
            file.write(f"{key} : {monitoring_values[key]}\n")

        file.write(f"compression rate : {compression_rate}")


if __name__ == "__main__":

    main()
