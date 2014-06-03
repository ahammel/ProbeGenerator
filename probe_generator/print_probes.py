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


class Nothing(object):
    """Represents a failed computation.

    """


class TryChain(object):
    """Apply a series of functions to a value, keeping only the first one that
    fails to throw a user-specified exception.

    E.g:

        >>> c = TryChain(ZeroDivisionError)
        >>> c.bind(lambda : 1 / '3')
        # Only ZeroDivisionError is caught. TypeError is thrown as normal.
        >>> c.bind(lambda : 1 / 0)
        >>> c.bind(lambda : 4 / 2)
        >>> c.bind(lambda : 5 / 3)
        >>> c.value
        2.0

    """
    def __init__(self, exception):
        self._exception = exception
        self.error = Nothing
        self.value = Nothing

    def bind(self, function):
        """If the value of the chain is Nothing, call the (nullary) function
        and apply it to the value.

        If the function call raises an exception, store the value of that
        exception in error property.

        """
        if self.value is Nothing:
            try:
                self.value = function()
            except self._exception as error:
                self.error = error


def print_probes(statement_file, genome_file, *annotation_files):
    """Print probes in FASTA format given a reference genome file and a file
    containing SNP probe statements.

    """
    with open(statement_file) as statements, open(genome_file) as genome:
        ref_genome = reference.reference_genome(genome)
        annotations = _combine_annotations(annotation_files)
        for statement in statements:
            chain = TryChain(InvalidStatement)
            chain.bind(lambda : [CoordinateProbe.from_statement(statement)])
            chain.bind(lambda : [SnpProbe.from_statement(statement)])
            chain.bind(lambda : ExonProbe.explode(statement, annotations))

            if chain.value is Nothing:
                raise chain.error

            probes = chain.value

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
