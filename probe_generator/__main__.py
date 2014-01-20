"""Automatically generate probe sequences.

Usage:
    probe-generator --coordinate COORD  --genome GENOME

Options:
    -c COORD --coordinate=COORD     a coordinate statement
                                    e.g. "1:100-10/2:200+20"
    -g GENOME --genome=GENOME       the Ensembl reference genome
                                    (FASTA format)

"""
from probe_generator import reference, coordinate_statement

from docopt import docopt

VERSION = '0.0'


def print_fasta(head, sequence):
    """Print a single string in FASTA format.

    """
    print(">{}\n{}\n".format(head, sequence))


def bases_from_coordinate(coordinate, ref_genome):
    """Given a set of coordinates and a genome, return a probe sequence.

    """
    first_bases = reference.bases(
            ref_genome,
            coordinate['chromosome1'],
            coordinate['start1'],
            coordinate['end1'])
    second_bases = reference.bases(
            ref_genome,
            coordinate['chromosome2'],
            coordinate['start2'],
            coordinate['end2'])
    return first_bases + second_bases


def print_probe_from_coordinate(statement, genome_file):
    """Print a probe in FASTA format given a coordinate statement and a
    reference genome file.

    """
    with open(genome_file) as handle:
        ref_genome = reference.reference_genome(handle)
        coordinate = coordinate_statement.parse(statement)
        bases = bases_from_coordinate(coordinate, ref_genome)
        print_fasta(statement, bases)


def main():
    args = docopt(__doc__, version='ProbeGenerator {}'.format(VERSION))
    print_probe_from_coordinate(args['--coordinate'], args['--genome'])


if __name__ == '__main__':
    main()
