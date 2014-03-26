"""Automatically generate probe sequences.

Usage:
    probe-generator --statement STMT  --genome GENOME --annotation FILE...  [-f]
    probe-generator --coordinate COORD  --genome GENOME [-f]

Options:
    -c COORD --coordinate=COORD     a file containing coordinate statements
                                    e.g. "1:100-10/2:200+20"
    -s STMT --statement=STMT        a file containing fusion statements
                                    e.g. "FOO#exon[1]-10/BAR#exon[2]+20"
    -g GENOME --genome=GENOME       the Ensembl reference genome
                                    (FASTA format)
    -a FILE --annotation=FILE       a genome annotation file in UCSC format
    -f --force                      run even if the total system memory is
                                    insufficient or cannot be determined

Please note that using --force is *strongly* discouraged.

"""
import sys

from docopt import docopt

from probe_generator import print_probes, check_memory

VERSION = '0.2.3'

REQUIRED_SYSTEM_MEMORY = 10485760 # 10Gb in Kb


def main():
    args = docopt(__doc__, version='ProbeGenerator {}'.format(VERSION))
    if not args['--force']:
        try:
            if check_memory.total_ram() < REQUIRED_SYSTEM_MEMORY:
                print("\nWARNING: total system memory is less than the "
                      "recommended minimum\n\n"
                      "Use '--force' to run anyway, but ONLY IF YOU KNOW "
                      "WHAT YOU'RE DOING\n\n"
                      "Alternatively, run probe-generator in in a "
                      "high-memory environment such as xhost08 or the\n"
                      "genesis cluster\n",
                      file=sys.stderr)
                sys.exit(1)
        except check_memory.Error as error:
            print("\nTotal system memory could not be determined: {}\n\n"
                  "Use '--force' to run anyway, but only if you have 16Gb "
                  "of RAM or so available, or you know what you're doing\n\n"
                  "See README.md for details".format(error),
                  file=sys.stderr)
            sys.exit(1)
    if args['--coordinate'] is not None:
        print_probes.from_coordinate(args['--coordinate'], args['--genome'])
    else:
        print_probes.from_statements(
                args['--statement'], args['--genome'], args['--annotation'])


if __name__ == '__main__':
    main()
