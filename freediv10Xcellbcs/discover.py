import os
import logging
import numpy as np
from collections import Counter
from .misc import load_bc_list
from .aligners import PrimerAligner 
from .plotting import knee_plot


log = logging.getLogger(__name__)


def discover_bcs_main(arguments):
    """
    Discover and demultiplex barcodes.
    """
    if not os.path.exists(arguments.output_dir):
        os.makedirs(arguments.output_dir)

    fwd_primer_max_end, rc_primer_max_end = find_primer_dist_ends(arguments)

    bc_oi_list = discover_bcs(fwd_primer_max_end, rc_primer_max_end, arguments)



def find_primer_dist_ends(arguments):
    # Find at least 10k seqs each with primer on fwd and rev strand
    # Return observed end based on percentiles with offset
    log.info('Finding primer ends...')
    fwd_ends, rc_ends = [], []
    primer_aligner = PrimerAligner()
    for total_seqs, (rec, strand, start, end) in enumerate(
            primer_aligner.iterate_recs_with_primer_pos(arguments.fastq_files)
            ):
        if strand == '+':
            fwd_ends.append(end)
        elif strand == '-':
            rc_ends.append(end)
        if len(fwd_ends) >= 10000 and len(rc_ends) >= 10000:
            break
    else:
        if len(fwd_ends) < 100 or len(rc_ends) < 100:
            raise RuntimeError(f'Not enough reads on fwd or rc strand: {len(fwd_ends):,d}, {len(rc_ends)}')
        log.warn(f'Did not find 10k fwd and rc seqs. Actual: {len(fwd_ends):,d}, {len(rc_ends)}')
    total_seqs += 1
    log.info(f'Found {len(fwd_ends):,d} fwd and {len(rc_ends):,d} rc primers after {total_seqs:,d} reads')
    fwd_90 = int(np.percentile(fwd_ends, 90))
    fwd_99 = int(np.percentile(fwd_ends, 99))
    rc_90 = int(np.percentile(rc_ends, 90))
    rc_99 = int(np.percentile(rc_ends, 99))
    fwd_end = min(fwd_90 + 20, fwd_99 + 8)
    rc_end = min(rc_90 + 20, rc_99 + 8)
    log.debug(f'Fwd ends 90th pctl = {fwd_90:,d}, 99th pctl = {fwd_99:,d}')
    log.debug(f'Rc ends 90th pctl = {rc_90:,d}, 99th pctl = {rc_99:,d}')
    log.info(f'Fwd primer search end: {fwd_end:,d},  Rc primer search end: {rc_end:,d}')
    return fwd_end, rc_end


def discover_bcs(fwd_primer_max_end, rc_primer_max_end, arguments):
    # Find primer ends and gather the raw 16bp following the end
    if arguments.barcode_whitelist:
        log.info('Loading whitelist...')
        whitelist = set(load_bc_list(arguments.barcode_whitelist))
        def in_whitelist(bc):
            return bc in whitelist
    else:
        def in_whitelist(bc):
            return True

    desired_bcs = 5000 * arguments.expected_cells
    bc_cntr = Counter()
    found_bcs = 0
    primer_aligner = PrimerAligner(fwd_primer_max_end, rc_primer_max_end)
    log.info('Discovering bcs...')
    for total_seqs, (rec, strand, start, end) in enumerate(
            primer_aligner.iterate_recs_with_primer_pos(arguments.fastq_files)
            ):
        if total_seqs % 1000000 == 0:
            log.info(f'\t{total_seqs:,d}')
        if strand is None:
            continue
        elif strand == '-':
            rec = rec.reverse_complement()
        bc = str(rec.seq)[end:end+16]

        if not in_whitelist(bc):
            continue
        bc_cntr[bc] += 1
        found_bcs += 1
        if found_bcs >= desired_bcs:
            break
    else:
        log.warn(f'Low barcode count: Only found {found_bcs:,d} of {desired_bcs:,d} desired raw bcs')

    bc_oi_thresh = 50
    bcs_and_counts = [(bc, count) for bc, count in bc_cntr.items()]
    bcs_and_counts.sort(reverse=True, key=lambda tup: tup[1])
    counts = [count for bc, count in bcs_and_counts]
    bc_oi_list = [bc for bc, count in bcs_and_counts if count > bc_oi_thresh]
    with open(os.path.join(arguments.output_dir, 'all_bcs_and_counts.txt'), 'w') as out:
        out.write('\n'.join([f'{bc}\t{count}' for bc, count in bcs_and_counts]))
    with open(os.path.join(arguments.output_dir, 'bcs_of_interest.txt'), 'w') as out:
        out.write('\n'.join(bc_oi_list))
    log.info(f'Found {len(bc_oi_list):,d} barcodes of interest')

    fig, ax = knee_plot(counts, bc_oi_thresh, good_label='BCs of interest')
    fig.savefig(os.path.join(arguments.output_dir, 'bcs_of_interest_knee_plot.png'), dpi=300)

    return bc_oi_list
