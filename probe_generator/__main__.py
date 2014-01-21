"""Automatically generate probe sequences.

Usage:
    probe-generator --statement STMT  --genome GENOME --annotation FILE...
    probe-generator --coordinate COORD  --genome GENOME

Options:
    -c COORD --coordinate=COORD     a file contatinng coordinate statements
                                    e.g. "1:100-10/2:200+20"
    -s STMT --statement=STMT        a file containg fusion statements
                                    e.g. "FOO#exon[1]-10/BAR#exon[2]+20"
    -g GENOME --genome=GENOME       the Ensembl reference genome
                                    (FASTA format)
    -a FILE --annotation=FILE       a genome annotation file in UCSC format

"""
import itertools

from probe_generator import (reference,
                             coordinate_statement,
                             annotation,
                             probe_statement,
                             sequence)

from docopt import docopt

VERSION = '0.0'

_COMPLEMENT = str.maketrans('acgtACGT', 'tgcaTGCA')


def reverse_complement(string):
    """Return the reverse-complement of a string of nucleotides.

    """
    return ''.join(reversed(string.translate(_COMPLEMENT)))


def print_fasta(head, bases):
    """Print a single string in FASTA format.

    """
    print(">{}\n{}".format(head, bases))


def combine_annotations(annotation_files):
    """Given a list of annotation files, return a single annotation.

    """
    rows = []
    for annotation_file in annotation_files:
        with open(annotation_file) as handle:
            rows.extend(annotation.parse_ucsc_file(handle))
    return rows


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
    if coordinate['inversion']:
        second_bases = reverse_complement(second_bases)
    return first_bases + second_bases


def print_probes_from_coordinate(statements_file, genome_file):
    """Print a probe in FASTA format given a reference genome file and a file
    containing coordinate statements.

    """
    with open(genome_file) as genome, open(statements_file) as statements:
        ref_genome = reference.reference_genome(genome)
        for statement in statements:
            statement = statement.strip()
            coordinate = coordinate_statement.parse(statement)
            bases = bases_from_coordinate(coordinate, ref_genome)
            print_fasta(statement, bases)


def print_probes_from_statements(statements_file, genome_file, annotation_files):
    """Print probes in FASTA format give a rerence genome file, a file
    containing probe statements, and UCSC genome annotations.

    """
    combined_annotation = combine_annotations(annotation_files)
    with open(statements_file) as statements, open(genome_file) as genome:
        ref_genome = reference.reference_genome(genome)
        for statement in statements:
            statement = statement.strip()
            specification = probe_statement.parse(statement)
            expanded_specifications = probe_statement.expand(specification)
            for expanded_specification in expanded_specifications:
                left_rows = annotation.lookup_gene(
                        specification['gene1'], combined_annotation)
                right_rows = annotation.lookup_gene(
                        specification['gene2'], combined_annotation)
                for left, right in itertools.product(left_rows, right_rows):
                    coordinate = sequence.sequence_range(
                            specification, left, right)
                    bases = bases_from_coordinate(coordinate, ref_genome)
                    header = "{} {} {}".format(
                            statement, right['name'], left['name'])
                    print_fasta(header, bases)


def main():
    args = docopt(__doc__, version='ProbeGenerator {}'.format(VERSION))
    if args['--coordinate'] is not None:
        print_probes_from_coordinate(args['--coordinate'], args['--genome'])
    else:
        print_probes_from_statements(
                args['--statement'], args['--genome'], args['--annotation'])


if __name__ == '__main__':
    main()
