import gzip

def gzip_friendly_open(fpath):
    return gzip.open(fpath, 'rt') if fpath.endswith('gz') else open(fpath)

def load_bc_list(fpath):
    bc_list = [line.strip() for line in gzip_friendly_open(fpath)]
    if not bc_list:
        raise RuntimeError(f'No barcodes found in {arguments.barcode_file}')
    bc_list = [bc[:bc.index('-')] if '-' in bc else bc for bc in bc_list] # clean cellranger bcs
    return bc_list
