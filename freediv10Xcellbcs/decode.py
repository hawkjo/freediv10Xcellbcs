import logging
import numpy as np
from Bio import SeqIO
from collections import Counter
from .misc import gzip_friendly_open, load_bc_list
from .aligners import PrimerAligner, OrientedPrimerSeq, BcUmiTailAligner
from .plotting import knee_plot


log = logging.getLogger(__name__)


def decode_fastqs(arguments):
    """
    Discover and demultiplex barcodes.
    """
    fwd_primer_max_end, rc_primer_max_end = find_primer_ends(arguments)

    if arguments.barcode_file is None:
        bc_oi_list = discover_bcs(fwd_primer_max_end, rc_primer_max_end, arguments)
    else:
        log.info('Loading given barcode file...')
        bc_oi_list = load_bc_list(arguments.barcode_file)

    demultiplex_bcs(fwd_primer_max_end, rc_primer_max_end, bc_oi_list, arguments)


def find_primer_ends(arguments):
    # Find at least 10k seqs each with primer on fwd and rev strand
    # Return observed end based on percentiles with offset
    log.info('Finding primer ends...')
    primer_aligner = PrimerAligner()
    err_thresh = -int(0.25*len(primer_aligner.primerseq))  # no more than 25% errors
    fwd_ends, rc_ends = [], []
    total_seqs = 0
    for fastq_file in arguments.fastq_files:
        for rec in SeqIO.parse(gzip_friendly_open(fastq_file)):
            total_seqs += 1
            best_fwd_alignment = primer_aligner(rec)[0]
            best_rc_alignment = primer_aligner(rec.reverse_complement())[0]
            if best_fwd_alignment.score >= err_thresh and best_rc_alignment.score < err_thresh:
                fwd_ends.append(best_fwd_alignment.aligned[1][-1][-1])  # end position on seq
            elif best_fwd_alignment.score < err_thresh and best_rc_alignment.score >= err_thresh:
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
    log.info('Discovering bcs...')
    primer_aligner = PrimerAligner()
    err_thresh = -int(0.25*len(primer_aligner.primerseq))  # no more than 25% errors
    bc_cntr = Counter()
    if arguments.barcode_whitelist:
        log.info('Loading whitelist...')
        whitelist = set(load_bc_list(arguments.barcode_whitelist))
        def in_whitelist(bc):
            return bc in whitelist
    else:
        def in_whitelist(bc):
            return True

    desired_bcs = 5000 * arguments.expected_cells
    total_seqs = 0
    found_bcs = 0
    for fastq_file in arguments.fastq_files:
        for rec in SeqIO.parse(gzip_friendly_open(fastq_file)):
            total_seqs += 1
            best_fwd_alignment = primer_aligner(rec[:fwd_primer_max_end])[0]
            best_rc_alignment = primer_aligner(rec.reverse_complement()[:rc_primer_max_end])[0]
            if best_fwd_alignment.score >= err_thresh and best_rc_alignment.score < err_thresh:
                fwd_end = best_fwd_alignment.aligned[1][-1][-1]
                bc = str(rec.seq)[fwd_end:fwd_end+16]
            elif best_fwd_alignment.score < err_thresh and best_rc_alignment.score >= err_thresh:
                rc_end = best_rc_alignment.aligned[1][-1][-1]
                bc = str(rec.reverse_complement().seq)[rc_end:rc_end+16]
            else:
                continue
            bc_cntr[bc] += 1
            found_bcs += 1
            if found_bcs >= desired_bcs:
                break
        else:
            continue
        break
    else:
        if found_bcs < desired_bcs / 2:
            raise RuntimeError(f'Less than half the desired amount of raw bcs found: {found_bcs:,d}')
        log.warn(f'Only found {found_bcs:,d} of {desired_bcs:,d} desired raw bcs')

    bc_oi_thresh = 50
    fig, ax = knee_plot(bc_cntr, bc_oi_thresh, good_label='BCs of interest')
    fig.savefig(os.path.join(arguments.output_dir, 'bcs_of_interest_knee_plot.png'), dpi=300)

    bcs_and_counts = [(bc, count) for bc, count in bc_cntr.items()]
    bcs_and_counts.sort(reverse=True, key=lambda tup: tup[1])
    bc_oi_list = [bc for bc, count in bcs_and_counts if count > bc_oi_thresh]
    with open(os.path.join(arguments.output_dir, 'all_bcs_and_counts.txt'), 'w') as out:
        out.write('\n'.join([f'{bc}\t{count}' for bc, count in bcs_and_counts]))
    with open(os.path.join(arguments.output_dir, 'bcs_of_interest.txt'), 'w') as out:
        out.write('\n'.join(bc_oi_list))
    log.info(f'Found {len(bc_oi_list):,d} barcodes of interest')

    return bc_oi_list
    

def demultiplex_bcs(fwd_primer_max_end, rc_primer_max_end, bc_oi_list, arguments):
    #TODO
