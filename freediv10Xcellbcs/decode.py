import logging
import numpy as np
from Bio import SeqIO
from collections import Counter
from .misc import gzip_friendly_open
from .aligners import PrimerAligner, OrientedPrimerSeq, BcUmiTailAligner


log = logging.getLogger(__name__)


def decode_fastqs(arguments):
    """
    Discover and demultiplex barcodes.
    """
    fwd_primer_max_end, rc_primer_max_end = find_primer_ends(arguments)

    if arguments.barcode_file is None:
        bc_oi_list = discover_bcs(fwd_primer_max_end, rc_primer_max_end, arguments)
    else:
        bc_oi_list = [line.strip() for line in gzip_friendly_open(arguments.barcode_file)]
        if not bc_oi_list:
            raise RuntimeError(f'No barcodes found in {arguments.barcode_file}')
        bc_oi_list = [bc[:bc.index('-')] if '-' in bc else bc for bc in bc_oi_list] # clean cellranger bcs

    demultiplex_bcs(fwd_primer_max_end, rc_primer_max_end, bc_oi_list, arguments)


def find_primer_ends(arguments):
    # Find at least 10k seqs each with primer on fwd and rev strand
    # Return observed end based on percentiles with offset
    log.info('Finding primer ends...')
    primer_aligner = PrimerAligner()
    thresh = -int(0.25*len(primer_aligner.primerseq))  # no more than 25% errors
    fwd_ends, rc_ends = [], []
    total_seqs = 0
    for fastq_file in arguments.fastq_files:
        for rec in SeqIO.parse(gzip_friendly_open(fastq_file)):
            total_seqs += 1
            best_fwd_alignment = primer_aligner(rec)[0]
            best_rc_alignment = primer_aligner(rec.reverse_complement())[0]
            if best_fwd_alignment.score >= thresh and best_rc_alignment.score < thresh:
                fwd_ends.append(best_fwd_alignment.aligned[1][-1][-1])  # end position on seq
            elif best_fwd_alignment.score < thresh and best_rc_alignment.score >= thresh:
                rc_ends.append(best_rc_alignment.aligned[1][-1][-1])
            if len(fwd_ends) >= 10000 and len(rc_ends) >= 10000:
                break
        else:
            continue
        break
    else:
        if len(fwd_ends) < 100 or len(rc_ends) < 100:
            raise RuntimeError(f'Not enough reads on fwd or rc strand: {len(fwd_ends):,d}, {len(rc_ends)}')
        log.warn(f'Did not find 10k fwd and rc seqs. Actual: {len(fwd_ends):,d}, {len(rc_ends)}')
    log.info(f'Found {len(fwd_ends):,d} fwd and {len(rc_ends)} rc primers after {total_seqs:,d} seqs')
    fwd_75 = np.percentile(fwd_ends, 75)
    fwd_99 = np.percentile(fwd_ends, 99)
    rc_75 = np.percentile(rc_ends, 75)
    rc_99 = np.percentile(rc_ends, 99)
    fwd_end = min(fwd_75 + 20, fwd_99 + 5)
    rc_end = min(rc_75 + 20, rc_99 + 5)
    log.debug(f'Fwd ends 75th pctl = {fwd_75:,d}, 99th pctl = {fwd_99:,d}')
    log.debug(f'Rc ends 75th pctl = {rc_75:,d}, 99th pctl = {rc_99:,d}')
    log.info(f'Fwd end: {fwd_end:,d},  Rc end: {rc_end:,d}')
    return fwd_end, rc_end


def discover_bcs(fwd_primer_max_end, rc_primer_max_end, arguments):
    # Find primer ends and gather the raw 16bp following the end
    # TODO

def demultiplex_bcs(fwd_primer_max_end, rc_primer_max_end, bc_oi_list, arguments):
    #TODO
