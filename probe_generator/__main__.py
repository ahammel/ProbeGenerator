"""Automatically generate probe sequences.

Usage:
    probe-generator --statements FILE --genome FILE [--annotation FILE...] [-f]

Options:
    -s FILE --statements=FILE       a file containing probe statements
    -g FILE --genome=FILE           the reference genome (FASTA format)
    -a FILE --annotation=FILE       a genome annotation file in UCSC format
    -f --force                      run even if the total system memory is
                                    insufficient or cannot be determined

"""
import sys

from docopt import docopt

from probe_generator import print_probes, check_memory

VERSION = '0.5'

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
    print_probes.print_probes(
            args['--statements'], args['--genome'], *args['--annotation'])


if __name__ == '__main__':
    main()
