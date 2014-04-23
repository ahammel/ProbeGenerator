"""Find the sequences of probes and print them in FASTA format.

"""
from probe_generator import reference, annotation
from probe_generator.coordinate_probe import CoordinateProbe
from probe_generator.snp_probe import SnpProbe
from probe_generator.exon_probe import ExonProbe
from probe_generator.probe import InvalidStatement


def print_probes(statement_file, genome_file, *annotation_files):
    """Print probes in FASTA format given a reference genome file and a file
    containing SNP probe statements.

    """
    with open(statement_file) as statements, open(genome_file) as genome:
        ref_genome = reference.reference_genome(genome)
        for statement in statements:
            try:
                probes = [CoordinateProbe.from_statement(statement)]
            except InvalidStatement:
                try:
                    probes = [SnpProbe.from_statement(statement)]
                except InvalidStatement:
                    annotations = _combine_annotations(annotation_files)
                    probes = ExonProbe.explode(statement, annotations)
            # TODO: would be nice to cook up a monad for the above hideousness.
            #
            # E.g.:
            #   probes =<< [CoordinateProbe.from_statement(statement)]
            #   probes =<< [SnpProbe.from_statement(statement)]
            #   probes =<< ExonProbe.explode(statement, annotations)
            for probe in probes:
                print_fasta(probe, probe.sequence(ref_genome))


def print_fasta(head, bases):
    """Print a single string in FASTA format.

    """
    print(">{}\n{}".format(head, bases))


def _combine_annotations(annotation_files):
    """Given a list of annotation files, return a single annotation.

    """
    rows = []
    for annotation_file in annotation_files:
        with open(annotation_file) as handle:
            rows.extend(annotation.parse_ucsc_file(handle))
    return rows
