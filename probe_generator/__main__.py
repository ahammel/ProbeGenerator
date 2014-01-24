"""Automatically generate probe sequences.

Usage:
    probe-generator --statement STMT  --genome GENOME --annotation FILE...
    probe-generator --coordinate COORD  --genome GENOME

Options:
    -c COORD --coordinate=COORD     a file containing coordinate statements
                                    e.g. "1:100-10/2:200+20"
    -s STMT --statement=STMT        a file containing fusion statements
                                    e.g. "FOO#exon[1]-10/BAR#exon[2]+20"
    -g GENOME --genome=GENOME       the Ensembl reference genome
                                    (FASTA format)
    -a FILE --annotation=FILE       a genome annotation file in UCSC format

"""
from docopt import docopt

from probe_generator import print_probes

VERSION = '0.0'


def main():
    args = docopt(__doc__, version='ProbeGenerator {}'.format(VERSION))
    if args['--coordinate'] is not None:
        print_probes.from_coordinate(args['--coordinate'], args['--genome'])
    else:
        print_probes.from_statements(
                args['--statement'], args['--genome'], args['--annotation'])


if __name__ == '__main__':
    main()
