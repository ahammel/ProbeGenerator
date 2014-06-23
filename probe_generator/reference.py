"""Parse and extract base pair sequences from an Ensembl reference genome.

"""
from probe_generator import sequence
from probe_generator.exceptions import NonFatalError

def bases(sequence_range, genome):
    """Return the bases from a SequenceRange object.

    """
    raw_bases = _raw_bases(
        sequence_range.chromosome,
        sequence_range.start,
        sequence_range.end,
        genome)
    if sequence_range.reverse_complement:
        return sequence.reverse_complement(raw_bases)
    else:
        return raw_bases


def reference_genome(genome):
    """Map chromosomes to base pair sequences.

    `genome` is a handle to a reference genome in Ensembl FASTA format.

    Returns a dictionary.

    """
    genome_map = {}
    chromosome = None
    for line in genome:
        if line.startswith('>'):
            chromosome = line[1:].split()[0]
            # In an Ensembl reference genome, the chromosome is the first
            # string of characters after the '>' but before whitespace.
            # E.g.:
            #   >chr Homo spaiens some chromosome etc etc
            #   NNN...
            genome_map[chromosome] = []
        elif chromosome is None:
            raise InvalidGenomeFile(
                    "could not parse input: {!r}".format(
                        line))
        else:
            genome_map[chromosome].append(line.strip())
    if not genome_map:
        raise InvalidGenomeFile("genome file empty!")
    return {chromosome: ''.join(bases)
            for (chromosome, bases)
            in genome_map.items()}


def _raw_bases(chromosome, start, end, genome):
    """Return a string of the base pairs of chromosome from start to end.

    The start and end attributes follow the Python convention for slices
    (indexed from zero, start inclusive, end exclusive).

    The genome is a dictionary relating chromosome names to base pair sequences
    (which are strings).

    """
    try:
        base_pairs = genome[chromosome][start:end]
    except KeyError:
        raise MissingChromosome(
                "no such chromosome: {!r}".format(
                    chromosome))
    if end - start != len(base_pairs):
        raise NonContainedRange(
                "range [{0}:{1}] outside the "
                "range of chromosome {2!r}".format(
                    start, end, chromosome))
    return base_pairs



class NonContainedRange(Exception):
    """Raised when the range of base pairs which is to be sliced from a
    chromosome includes base pairs outside the chromosome.

    """


class InvalidGenomeFile(Exception):
    """Raised when a a genome_file cannot be parsed.

    """


class MissingChromosome(NonFatalError):
    """Raised when a chromosome is not present in the reference genome.

    """
