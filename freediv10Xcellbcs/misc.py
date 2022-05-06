import gzip

def gzip_friendly_open(fpath):
    return gzip.open(fpath, 'rt') if fpath.endswith('gz') else open(fpath)

