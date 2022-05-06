"""
Free divergence-based decoding of 10X cell barcodes

Usage:
  freediv10Xcellbcs decode       <fastq_files> <5p_or_3p> [--barcode-file=<barcode_file>] [--barcode-whitelist=<barcode_whitelist>] [--output-dir=<output_dir>] [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  decode        Decode barcodes in fastq files with same barcodes. Separate file names with commas.

"""
import logging
import os
from docopt import docopt
from .__init__ import __version__
from .config import CommandLineArguments
from .decode import decode_fastqs


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
    }

    commands[arguments.command](arguments)


if __name__ == '__main__':
    main()
