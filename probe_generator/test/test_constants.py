"""Constants shared among test modules.

"""
import os

VALIDATION_DATA_DIR = os.path.join(
        "/",
        "genesis",
        "scratch",
        "validations",
        "test_probe_generator")

GENOME =  {
        "1": "acgtacgt",
        "2": "aaaagggg",
        }

ANNOTATION = [
    {'name':       'FOO',
     'proteinID':  'ABC',
     'exonStarts': '1,3,',
     'exonEnds':   '2,4,',
     'strand':     '+',
     'chrom':      '1',
     },
    {'name':       'BAR',
     'proteinID':  'DEF',
     'exonStarts': '7,9,11,',
     'exonEnds':   '8,10,12,',
     'strand':     '+',
     'chrom':      '2',
     },
    {'name':       'BAZ',
     'proteinID':  'GHI',
     'exonStarts': '10,20',
     'exonEnds':   '15,25',
     'strand':     '-',
     'chrom':      '3',
     },
    ]
