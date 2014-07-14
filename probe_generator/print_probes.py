"""Find the sequences of probes and print them in FASTA format.

"""
import sys

# Utilities
from probe_generator import reference, annotation
# Probe classes
from probe_generator.coordinate_probe import CoordinateProbe
from probe_generator.snp_probe        import SnpProbe
from probe_generator.gene_snp_probe   import GeneSnpProbe
from probe_generator.amino_acid_probe import AminoAcidProbe
from probe_generator.exon_probe       import ExonProbe
# Exceptions
from probe_generator.probe import InvalidStatement, NonFatalError

NO_PROBES_WARNING = (
    "WARNING: no probes could be generated for statement {!r}")

INVALID_STATEMENT_WARNING = (
    "WARNING: the statement {!r} could not be parsed")


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

    def bind_all(self, *functions):
        """Apply `bind` to all the arguments in order.

        """
        for function in functions:
            self.bind(function)


def print_probes(statement_file, genome_file, *annotation_files):
    """Print probes in FASTA format given a reference genome file and a file
    containing SNP probe statements.

    """
    with open(statement_file) as statements, open(genome_file) as genome:
        ref_genome = reference.reference_genome(genome)
        annotations = _combine_annotations(annotation_files)
        for statement in statements:
            chain = TryChain(InvalidStatement)
            chain.bind_all(
                lambda: list(CoordinateProbe.explode(statement)),
                lambda: list(SnpProbe.explode(statement)),
                lambda: list(GeneSnpProbe.explode(statement, annotations)),
                lambda: list(AminoAcidProbe.explode(statement, annotations)),
                lambda: list(ExonProbe.explode(statement, annotations)))

            if chain.value is Nothing:
                print(INVALID_STATEMENT_WARNING.format(statement),
                      file=sys.stderr)
            probes = chain.value

            one_probe_printed = False
            for probe in probes:
                try:
                    print_fasta(probe, probe.sequence(ref_genome))
                except NonFatalError as error:
                    print(
                        "In probe: {}: {}".format(probe, error),
                        file=sys.stderr)
                else:
                    one_probe_printed = True

            if not one_probe_printed: # i.e., the generator was empty
                print(NO_PROBES_WARNING.format(statement), file=sys.stderr)


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
