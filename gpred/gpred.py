import argparse
import sys
import os
import csv
import re
import textwrap
from re import Pattern
from pathlib import Path
from typing import List, Union, Optional
__author__ = "Desvilles Aurélien"
__copyright__ = "Universite Paris Cité"
__credits__ = ["Desvilles Aurélien"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Desvilles Aurélien"
__email__ = "desvillesaurelien@gmail.com"
__status__ = "Done"


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


def get_arguments(): # pragma: no cover
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
    try:
        with open(fasta_file, "r") as file:
            lines = file.readlines()
            for line in lines:
                line = line.strip()  
                if line.startswith(">"):
                    continue  
                else:
                    sequence += line.upper()  
        return sequence
    except FileNotFoundError:
        print(f"Le fichier '{fasta_file}' est introuvable.")
    except Exception as e:
        print(f"Une erreur s'est produite : {e}")


def find_start(start_regex: Pattern, sequence: str, start: int, stop: int) -> Union[int, None]:
    """Find next start codon before a end position.

    :param start_regexp: A regex object that identifies a start codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :param stop: (int) Stop position of the research
    :return: (int) If exist, position of the start codon. Otherwise None. 
    """
    
    start_regex = re.compile(start_regex)

    match = start_regex.search(sequence, start, stop)
    if match:
        return match.start(0)
    else:
        return None


def find_stop(stop_regex: Pattern, sequence: str, start: int) -> Union[int, None]:
    """Find next stop codon that should be in the same reading phase as the start.

    :param stop_regexp: A regex object that identifies a stop codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :return: (int) If exist, position of the stop codon. Otherwise None. 
    """
    
    stop_regex = re.compile(stop_regex)

    matches = stop_regex.finditer(sequence, start)

    is_match = None
    for match in matches:
        if (match.end(0) - start) % 3 == 0:
            is_match = match.start(0)
            break
        else:
            is_match = None
    return is_match

def has_shine_dalgarno(shine_regex: Pattern, sequence: str, start: int, max_shine_dalgarno_distance: int) -> bool:
    """Find a shine dalgarno motif before the start codon

    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Position of the start in the genome
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno to the start position
    :return: (boolean) true -> has a shine dalgarno upstream to the gene, false -> no
    """
    
    shine_dalgarno_regex = re.compile(shine_regex)
    
    is_sd = False
    if start - max_shine_dalgarno_distance > 0:
        if (start - 6) > (start - max_shine_dalgarno_distance):
            match = shine_dalgarno_regex.search(sequence, start - max_shine_dalgarno_distance, start - 6)
   
            if match :
                is_sd = True
            else:
                is_sd = False
    
    return is_sd




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
    genes = []
    position_courante = 0
    longueur_sequence = len(sequence)

    while longueur_sequence - position_courante >= min_gap:
        start_position = find_start(start_regex, sequence, position_courante, longueur_sequence)
        
        if start_position is not None:
            stop_position = find_stop(stop_regex, sequence, start_position)
            
            if stop_position is not None:
                gene_length = stop_position - start_position 

                if gene_length >= min_gene_len:
                    shine_position = find_start(shine_regex, sequence, start_position - max_shine_dalgarno_distance, start_position)

                    if shine_position is not None and shine_position < start_position:
                        genes.append([start_position + 1, stop_position + 3])  

                        position_courante = stop_position + min_gap + 3
                    else:
                        position_courante = start_position + 1     
                else:
                    position_courante = start_position + 1       
            else:
                position_courante = start_position + 1
        else:
            position_courante = position_courante + 1

    return genes



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
                    (sequence[gene_pos[0]-1:gene_pos[1]])))
            i = i+1
            for j,gene_pos in enumerate(probable_genes_comp):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                            i+1+j, os.linesep,
                            (sequence_rc[gene_pos[0]-1:gene_pos[1]])))
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
def main() -> None: # pragma: no cover
    """
    Main program function
    """
    # Gene detection over genome involves to consider a thymine instead of
    # an uracile that we would find on the expressed RNA
    #start_codons = ['TTG', 'CTG', 'ATT', 'ATG', 'GTG']
    #stop_codons = ['TAA', 'TAG', 'TGA']
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    # Shine AGGAGGUAA
    #AGGA ou GGAGG 
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    # Arguments
    args = get_arguments()
    # Let us do magic in 5' to 3'
    
    # Don't forget to uncomment !!!
    # Call these function in the order that you want
    # We reverse and complement
    #sequence_rc = reverse_complement(sequence)
    # Call to output functions
    #write_genes_pos(args.predicted_genes_file, probable_genes)
    #write_genes(args.fasta_file, sequence, probable_genes, sequence_rc, probable_genes_comp)

    genome_sequence = read_fasta(args.genome_file)

    # Genes predicted in the 5' to 3' direction
    predicted_genes_5to3 = predict_genes(
        genome_sequence,
        start_regex,
        stop_regex,
        shine_regex,
        args.min_gene_len,
        args.max_shine_dalgarno_distance,
        args.min_gap
    )

    # Complementary sequence in direction 3' to 5'
    genome_sequence_rc = reverse_complement(genome_sequence)

    # Genes predicted in the 3' to 5' direction
    predicted_genes_3to5 = predict_genes(
        genome_sequence_rc,
        start_regex,
        stop_regex,
        shine_regex,
        args.min_gene_len,
        args.max_shine_dalgarno_distance,
        args.min_gap
    )

    # Correction of gene positions in the 5' to 3' direction
    for gene in predicted_genes_3to5:
        gene[0], gene[1] = len(genome_sequence) - gene[1], len(genome_sequence) - gene[0]

    # Merging results into a single list
    all_predicted_genes = predicted_genes_5to3 + predicted_genes_3to5

    # Sorting genes by start position
    all_predicted_genes.sort(key=lambda x: x[0])

    # Writing predicted gene positions
    write_genes_pos(args.predicted_genes_file, all_predicted_genes)

    # Writing predicted gene sequences
    write_genes(args.fasta_file, genome_sequence, all_predicted_genes, genome_sequence_rc, all_predicted_genes)

if __name__ == '__main__':
    main()
