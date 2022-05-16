"""
Free divergence-based decoding of 10X cell barcodes

Usage:
  freediv10Xcellbcs decode       <fastq_files> <kit_5p_or_3p> [--barcode-file=<barcode_file>] [--max-err-decode=<max_err_decode>] [--reject-delta=<reject_delta>] [--decoder-file=<decoder_file>] [--output-dir=<output_dir>] [--no-umis] [-v | -vv | -vvv]
  freediv10Xcellbcs discover     <fastq_files> [--barcode-whitelist=<barcode_whitelist>] [--expected-cells=<expected_cells>] [--output-dir=<output_dir>] [-v | -vv | -vvv]
  freediv10Xcellbcs prebuild     --barcode-file=<barcode_file> --max-err-decode=<max_err_decode> --reject-delta=<reject_delta> [--output-dir=<output_dir>] [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  decode        Decode barcodes in fastq files with same barcodes. Separate file names with commas.
  discover      Discover barcodes de novo in fastq files. Separate file names with commas.
  prebuild      Prebuild decoder for repeated use with same barcodes and settings.

"""
import logging
import os
from docopt import docopt
from .__init__ import __version__
from .config import CommandLineArguments
from .decode import decode_fastqs, build_decoder 
from .discover import discover_bcs_main


def main(**kwargs):
    docopt_args = docopt(__doc__, version=__version__)
    arguments = CommandLineArguments(docopt_args, os.getcwd())

    log = logging.getLogger()
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s   %(message)s", "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(arguments.log_level)
    log.debug(docopt_args)

    commands = {
        'decode': decode_fastqs,
        'discover': discover_bcs_main,
        'prebuild': build_decoder, 
    }

    commands[arguments.command](arguments)


if __name__ == '__main__':
    main()
