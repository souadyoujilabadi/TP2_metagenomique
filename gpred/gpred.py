import argparse
import sys
import os
import csv
import re
import textwrap
from re import Pattern
from pathlib import Path
from typing import List, Union, Optional


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments():  # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='genome_file', type=isfile, required=True, 
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', dest='min_gene_len', type=int, 
                        default=50, help="Minimum gene length to consider (default 50).")
    parser.add_argument('-s', dest='max_shine_dalgarno_distance', type=int, 
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif (default 16).")
    parser.add_argument('-d', dest='min_gap', type=int, default=40,
                        help="Minimum gap between two genes - shine box not included (default 40).")
    parser.add_argument('-p', dest='predicted_genes_file', type=Path, 
                        default=Path("predict_genes.csv"),
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=Path,
                        default=Path("genes.fna"),
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file: Path) -> str:
    """Extract genome sequence from fasta files.

    :param fasta_file: (Path) Path to the fasta file.
    :return: (str) Sequence from the genome. 
    """
    sequence = ""

    with fasta_file.open('r') as f:
        for line in f:
            # Skip the header
            if not line.startswith('>'):
                # Remove any whitespace and append to the sequence
                sequence += line.strip()

    return sequence.upper()


def find_start(start_regex: Pattern, sequence: str, start: int, stop: int) -> Union[int, None]:
    """Find next start codon before a end position.

    :param start_regexp: A regex object that identifies a start codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :param stop: (int) Stop position of the research
    :return: (int) If exist, position of the start codon. Otherwise None.
    """
    # Use the search method of the regex object to find the first occurrence
    # of a start codon between the given start and stop positions
    match = start_regex.search(sequence, start, stop)
    # If a match was found, return the start position of the match
    if match:
        return match.start()
    # If no match was found, return None
    return None


def find_stop(stop_regex: Pattern, sequence: str, start: int) -> Union[int, None]:
    """Find next stop codon that should be in the same reading phase as the start.

    :param stop_regexp: A regex object that identifies a stop codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :return: (int) If exist, position of the stop codon. Otherwise None. 
    """
    # Check every third position to maintain reading frame
    for i in range(start, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        # Check for stop codon match
        if stop_regex.match(codon):
            return i
    # Return None if no stop codon found
    return None


def has_shine_dalgarno(shine_regex: Pattern, sequence: str, start: int, max_shine_dalgarno_distance: int) -> bool:
    """Find a shine dalgarno motif before the start codon

    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Position of the start in the genome
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno to the start position
    :return: (boolean) true -> has a shine dalgarno upstream to the gene, false -> no
    """
    # Calculate the start and end positions for the search
    start_pos = start - max_shine_dalgarno_distance
    end_pos = start - 6
    return shine_regex.search(sequence, start_pos, end_pos) is not None if start_pos >= 0 else False


def predict_genes(sequence: str, start_regex: Pattern, stop_regex: Pattern, shine_regex: Pattern, 
                  min_gene_len: int, max_shine_dalgarno_distance: int, min_gap: int) -> List[List[int]]:
    """Predict most probable genes

    :param sequence: (str) Sequence from the genome.
    :param start_regexp: A regex object that identifies a start codon.
    :param stop_regexp: A regex object that identifies a stop codon.
    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param min_gene_len: (int) Minimum gene length.
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno to the start position.
    :param min_gap: (int) Minimum distance between two genes.
    :return: (list) List of [start, stop] position of each predicted genes.
    """
    # Initialize the current position and the list of predicted genes
    current_position = 0
    predicted_genes = []

    while len(sequence) - current_position >= min_gap:
        # Using find_start function to find the next start codon from current_position
        start = find_start(start_regex, sequence, current_position, len(sequence))
        if start is not None:
            # Using find_stop function to find the next stop codon from start position
            stop = find_stop(stop_regex, sequence, start)
            if stop is not None:
                # Check if the identified gene meets the minimum length criterion
                if stop - start + 3 >= min_gene_len:  # +3 to include the stop codon
                    # Check for Shine-Dalgarno sequence upstream of the start codon using has_shine_dalgarno function
                    if has_shine_dalgarno(shine_regex, sequence, start, max_shine_dalgarno_distance):
                        # Add the identified gene to the list of predicted genes
                        predicted_genes.append([start + 1, stop + 3])  # +1 because genomic position starts at 1; +3 to include the stop codon
                        # Update the current position to after the found stop codon
                        current_position = stop + 3 + min_gap
                        continue  # move to the next iteration
            current_position += 1
        else:
            current_position += 1

    return predicted_genes


def write_genes_pos(predicted_genes_file: Path, probable_genes: List[List[int]]) -> None:
    """Write list of gene positions.

    :param predicted_genes_file: (Path) Output file of gene positions.
    :param probable_genes: List of [start, stop] position of each predicted genes.
    """
    try:
        with predicted_genes_file.open("wt") as predict_genes:
            predict_genes_writer = csv.writer(predict_genes, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
    except IOError:
        sys.exit("Error cannot open {}".format(predicted_genes_file))


def write_genes(fasta_file: Path, sequence: str, probable_genes: List[List[int]], sequence_rc: str, 
                probable_genes_comp: List[List[int]]):
    """Write gene sequence in fasta format

    :param fasta_file: (Path) Output fasta file.
    :param sequence: (str) Sequence of genome file in 5'->3'.
    :param probable_genes: (list) List of [start, stop] position of each predicted genes in 5'->3'.
    :param sequence_rc: (str) Sequence of genome file in 3' -> 5'.
    :param probable_genes_comp: (list)List of [start, stop] position of each predicted genes in 3' -> 5'.
    """
    try:
        with open(fasta_file, "wt") as fasta:
            for i,gene_pos in enumerate(probable_genes):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                    i+1, os.linesep, 
                    textwrap.fill(sequence[gene_pos[0]-1:gene_pos[1]])))
            i = i+1
            for j,gene_pos in enumerate(probable_genes_comp):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                            i+1+j, os.linesep,
                            textwrap.fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])))
    except IOError:
        sys.exit("Error cannot open {}".format(fasta_file))


def reverse_complement(sequence: str) -> str:
    """Get the reverse complement

    :param sequence: (str) DNA Sequence.
    :return: (str) Reverse complemented sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in sequence[::-1]])


#==============================================================
# Main program
#==============================================================
def main() -> None:  # pragma: no cover
    """
    Main program function
    """
    # Gene detection over genome involves to consider a thymine instead of
    # an uracile that we would find on the expressed RNA
    # start_codons = ['TTG', 'CTG', 'ATT', 'ATG', 'GTG']
    # stop_codons = ['TAA', 'TAG', 'TGA']
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    # Shine AGGAGGUAA
    # AGGA ou GGAGG
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    # Arguments
    args = get_arguments()
    # Let us do magic in 5' to 3'

    # Read the genome sequence from the fasta file
    sequence = read_fasta(args.genome_file)

    # Predict probable genes in 5' to 3' direction
    probable_genes = predict_genes(sequence, start_regex, stop_regex, shine_regex,
                                   args.min_gene_len, args.max_shine_dalgarno_distance, args.min_gap)

    # Get the reverse complement of the sequence for 3' to 5' direction
    sequence_rc = reverse_complement(sequence)

    # Predict probable genes in 3' to 5' direction
    probable_genes_comp = predict_genes(sequence_rc, start_regex, stop_regex, shine_regex,
                                        args.min_gene_len, args.max_shine_dalgarno_distance, args.min_gap)

    # Correct the positions for genes predicted in 3' to 5' direction to appear in 5' to 3' orientation
    seq_len = len(sequence)
    probable_genes_comp_corrected = [[seq_len - pos[1] + 1, seq_len - pos[0] + 1] for pos in probable_genes_comp]

    # Merge the predicted genes from both directions
    all_predicted_genes = probable_genes + probable_genes_comp_corrected

    # Write the predicted gene positions and sequences to respective files
    write_genes_pos(args.predicted_genes_file, all_predicted_genes)
    write_genes(args.fasta_file, sequence, probable_genes, sequence_rc, probable_genes_comp_corrected)


if __name__ == '__main__':
    main()
