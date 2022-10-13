# Free divergence-based decoding of 10X cell barcodes

A package for decoding 10X barcodes with conflict-aware and rejection-shell capable FREE divergence DNA barcode error-correction as described in the manuscript:

### Manuscript to appear soon...

### Installation

The following instructions should work across platforms, except that installing virtualenv with apt-get is Ubuntu specific. For other platforms, install virtualenv appropriately if desired.

First, clone the repository to a local directory:

```
git clone https://github.com/hawkjo/freediv10Xcellbcs.git
```

Optionally, you can install into a virtual environment (recommended):

```
sudo apt-get install -y virtualenv
cd freediv10Xcellbcs
virtualenv envfreebarcodes
. envfreebarcodes/bin/activate
```

Now install required packages and freediv10Xcellbcs using pip:

```
pip install numpy==1.13.3 && pip install -r requirements.txt && python setup.py install
```

### Usage

```
Usage:
  freediv10Xcellbcs decode       <fastq_files> <kit_5p_or_3p> [--barcode-file=<barcode_file>] [--max-err-decode=<max_err_decode>] [--reject-delta=<reject_delta>] [--decoder-file=<decoder_file>] [--output-dir=<output_dir>] [--no-umis] [-v | -vv | -vvv]
  freediv10Xcellbcs discover     <fastq_files> [--barcode-whitelist=<barcode_whitelist>] [--expected-cells=<expected_cells>] [--reads-per-cell=<reads_per_cell>] [--threshold=<threshold>] [--output-dir=<output_dir>] [-v | -vv | -vvv]
  freediv10Xcellbcs prebuild     --barcode-file=<barcode_file> --max-err-decode=<max_err_decode> --reject-delta=<reject_delta> [--output-dir=<output_dir>] [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  decode        Decode barcodes in fastq files with same barcodes. Separate file names with commas.
  discover      Discover barcodes de novo in fastq files. Separate file names with commas.
  prebuild      Prebuild decoder for repeated use with same barcodes and settings.
```


