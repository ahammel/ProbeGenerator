"""Find the sequences of probes and print them in FASTA format.

"""
import sys

from probe_generator import reference, annotation
from probe_generator.coordinate_probe import CoordinateProbe
from probe_generator.snp_probe import SnpProbe
from probe_generator.exon_probe import ExonProbe
from probe_generator.probe import InvalidStatement, NonFatalError

NO_PROBES_WARNING = (
    "WARNING: no probes could be generated for statement {}\n\n"
    "This is usually beacuse exon probes are in use and no genome annotation"
    "file was provided\n\n"
    "See `probe-generator --help")


def print_probes(statement_file, genome_file, *annotation_files):
    """Print probes in FASTA format given a reference genome file and a file
    containing SNP probe statements.

    """
    with open(statement_file) as statements, open(genome_file) as genome:
        ref_genome = reference.reference_genome(genome)
        # See the 'monad' branch for a planned refactoring of this
        # try/excpet tree --- AJH
        for statement in statements:
            try:
                probes = [CoordinateProbe.from_statement(statement)]
            except InvalidStatement:
                try:
                    probes = SnpProbe.explode(statement)
                except InvalidStatement:
                    annotations = _combine_annotations(annotation_files)
                    probes = ExonProbe.explode(statement, annotations)

            if not probes:
                print(NO_PROBES_WARNING.format(statement), file=sys.stderr)

            for probe in probes:
                try:
                    print_fasta(probe, probe.sequence(ref_genome))
                except NonFatalError as error:
                    print(error, file=sys.stderr)


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

