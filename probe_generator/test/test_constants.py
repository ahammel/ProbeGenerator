"""Constants shared among test modules.

"""
import os

from probe_generator.transcript import Transcript

VALIDATION_DATA_DIR = os.path.join(
        "/",
        "genesis",
        "scratch",
        "validations",
        "test_probe_generator")

GENOME =  {
        "1": "acgtacgt",
        "2": "aaaaggggacata",
        "3": "aaaaaaaaacccGGGcccaaaaaaaaaaaaaaa",
        }

ANNOTATION = [
    Transcript(row) for row in (
        {'name':       'FOO',
         'proteinID':  'ABC',
         'exonStarts': '1,3,',
         'exonEnds':   '2,4,',
         'strand':     '+',
         'chrom':      '1',
         'cdsStart':   '1',
         'cdsEnd':     '4',
         },
        {'name':       'BAR',
         'proteinID':  'DEF',
         'exonStarts': '7,9,11,',
         'exonEnds':   '8,10,12,',
         'strand':     '+',
         'chrom':      '2',
         'cdsStart':   '7',
         'cdsEnd':     '12',
         },
        {'name':       'BAZ',
         'proteinID':  'GHI',
         'exonStarts': '10,20',
         'exonEnds':   '15,25',
         'strand':     '-',
         'chrom':      '3',
         'cdsStart':   '10',
         'cdsEnd':     '23',
         },
        {'name':       'BAZ2',
         'proteinID':  'GHI2',
         'exonStarts': '9',
         'exonEnds':   '18,',
         'strand':     '+',
         'chrom':      '3',
         'cdsStart':   '9',
         'cdsEnd':     '18',
         },
        {'name':       'QUX',
         'proteinID':  'JKL',
         'exonStarts': '100,200',
         'exonEnds':   '150,250',
         'strand':     '+',
         'chrom':      '4',
         'cdsStart':   '120',
         'cdsEnd':     '240',
         },
        {'name':       'FROB',
         'proteinID':  'MNO',
         'exonStarts': '6,12,18,',
         'exonEnds':   '9,15,21,',
         'strand':     '+',
         'chrom':      '3',
         'cdsStart':   '6',
         'cdsEnd':     '21',
         },
        )
    ]
