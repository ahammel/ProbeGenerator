"""Parse and extract base pair sequences from an Ensembl reference genome.

"""
from probe_generator import sequence, probe


def bases(ref_genome, chromosome, start, end):
    """Return the base pairs of `chromosome` from `start` to `end`.

    The string of base pairs is from the `start` to `end` _inclusive_,
    where the first base of the chromosome is 1.

    `ref_genome` is a dictionary relating chromosome names to base pair
    sequences (which are strings). `chromosome` is a string, `start` and `end`
    are ints.

    Returns a string of base pairs.

    """
    try:
        base_pairs = ref_genome[chromosome][start-1:end]
    except KeyError:
        raise MissingChromosome(
                "no such chromosome: {!r}".format(
                    chromosome))
    except TypeError:
        raise InvalidRange(
                "unsupported values for `start` and `end`: "
                "({!r}, {!r})".format(
                    start, end))
    if start > end:
        raise InvalidRange("unsupported values for `start` and `end`: "
                           "({}, {}). `start` must be  <= `end`".format(
                               start, end))
    if end-start+1 != len(base_pairs):
        raise NonContainedRange(
                "range [{0}:{1}] outside the "
                "range of chromosome {2!r}".format(
                    start, end, chromosome))
    return base_pairs


def bases_from_coordinate(coordinate, ref_genome):
    """Given a set of coordinates and a genome, return a probe sequence.

    """
    first_bases = bases(
            ref_genome,
            coordinate['chromosome1'],
            coordinate['start1'],
            coordinate['end1'])
    second_bases = bases(
            ref_genome,
            coordinate['chromosome2'],
            coordinate['start2'],
            coordinate['end2'])
    if coordinate['rc_side_1']:
        first_bases = sequence.reverse_complement(first_bases)
    if coordinate['rc_side_2']:
        second_bases = sequence.reverse_complement(second_bases)
    return first_bases + second_bases


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


class MissingChromosome(probe.NonFatalError):
    """Raised when a chromosome is not present in the reference genome.

    """


class InvalidRange(TypeError):
    """Raised on a TypeError while slicing a genome.

    """


class NonContainedRange(Exception):
    """Raised when the range of base pairs which is to be sliced from a
    chromosome includes base pairs outside the chromosome.

    """


class InvalidGenomeFile(Exception):
    """Raised when a a genome_file cannot be parsed.

    """
