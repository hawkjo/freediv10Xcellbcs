import os
import logging
import numpy as np
import freebarcodes
import freebarcodes.decode
from Bio import SeqIO
from collections import Counter
from .misc import load_bc_list, gzip_friendly_open, write_stats_file_from_cntr
from .aligners import PrimerAligner, OrientedPrimerSeq, BcUmiTailAligner
from .plotting import knee_plot, umi_len_plot


log = logging.getLogger(__name__)


def decode_fastqs(arguments):
    """
    Discover and demultiplex barcodes.
    """
    if not os.path.exists(arguments.output_dir):
        os.makedirs(arguments.output_dir)

    fwd_primer_max_end, rc_primer_max_end = find_primer_dist_ends(arguments)

    if arguments.barcode_file is None:
        bc_oi_list = discover_bcs(fwd_primer_max_end, rc_primer_max_end, arguments)
    else:
        log.info('Loading given barcode file...')
        bc_oi_list = load_bc_list(arguments.barcode_file)

    demultiplex_bcs(fwd_primer_max_end, rc_primer_max_end, bc_oi_list, arguments)


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
    fwd_75 = int(np.percentile(fwd_ends, 75))
    fwd_99 = int(np.percentile(fwd_ends, 99))
    rc_75 = int(np.percentile(rc_ends, 75))
    rc_99 = int(np.percentile(rc_ends, 99))
    fwd_end = min(fwd_75 + 20, fwd_99 + 8)
    rc_end = min(rc_75 + 20, rc_99 + 8)
    log.debug(f'Fwd ends 75th pctl = {fwd_75:,d}, 99th pctl = {fwd_99:,d}')
    log.debug(f'Rc ends 75th pctl = {rc_75:,d}, 99th pctl = {rc_99:,d}')
    log.info(f'Fwd primer search end: {fwd_end:,d},  Rc primer search end: {rc_end:,d}')
    return fwd_end, rc_end


def discover_bcs(fwd_primer_max_end, rc_primer_max_end, arguments):
    # Find primer ends and gather the raw 16bp following the end
    log.info('Discovering bcs...')
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
    for total_seqs, (rec, strand, start, end) in enumerate(
            primer_aligner.iterate_recs_with_primer_pos(arguments.fastq_files)
            ):
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
        if found_bcs < desired_bcs / 2:
            raise RuntimeError(f'Less than half the desired amount of raw bcs found for discovery: {found_bcs:,d} of {desired_bcs:,d}')
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
    

def demultiplex_bcs_and_umis(fwd_primer_max_end, rc_primer_max_end, bc_oi_list, arguments):
    def make_base_out_fpath(fpath):
        fname = os.path.basename(fpath)
        for ending in ['.fastq', '.fq', '.fastq.gz', '.fq.gz']:
            if fname.endswith(ending):
                fname = fname[:-len(ending)] + '.oriented_and_demult.d{}+{}'.format(
                        arguments.max_err_decode,
                        arguments.reject_delta
                        )
                return os.path.join(arguments.output_dir, fname)
        raise ValueError(f'Unrecognized fastq ending in {fpath}')

    log.info('Loading barcode decoder...')
    bd = freebarcodes.decode.FreeDivBarcodeDecoder()
    bd.build_codebook_from_random_codewords(bc_oi_list, arguments.max_err_decode, arguments.reject_delta)

    log.info('Building barcode-specific aligners...')
    bc_umi_aligners = {bc: BcUmiTailAligner(bc, 10, arguments.kit_5p_or_3p) for bc in bc_oi_list}
    example_full_prefix = bc_umi_aligners[bc_oi_list[0]].full_prefix

    log.info('Demultiplexing cell barcodes...')
    cum_stats = Counter()
    cum_umi_len_cntr = Counter()
    for fastq_fpath in arguments.fastq_files:
        log.info(f'Demultiplexing {fastq_fpath}')
        base_out_fpath = make_base_out_fpath(fastq_fpath)
        out_fastq_fpath = base_out_fpath + '.fastq'
        fq_out = open(out_fastq_fpath, 'w')
        log.info(f'Writing to {out_fastq_fpath}')
        stats = Counter()
        umi_len_cntr = Counter()
        for i, rec in enumerate(SeqIO.parse(gzip_friendly_open(fastq_fpath), 'fastq')):
            if i % 1000000 == 0:
                log.debug(f'\t{i:,d}')

            ops = OrientedPrimerSeq(rec, fwd_primer_max_end, rc_primer_max_end)
            if ops.orientation_failed:
                stats['Filter 1: orientation failed'] += 1
                continue
            stats['Filter 1: orientation successful'] += 1
            if not ops.is_well_formed:
                stats['Filter 2: primer alignment fails check'] += 1
                continue
            stats['Filter 2: primer alignment passes check'] += 1
            if len(ops.seq_after_primer) < len(example_full_prefix):
                stats['Filter 3: seq too short'] += 1
                continue
            stats['Filter 3: seq long enough'] += 1
            for obs_bc, start_pos in ops.get_kmers_at_primer_end_plusminus(16, 2):
                bc = bd.decode(obs_bc)
                if isinstance(bc, str):
                    ops.primer_end = start_pos
                    stats['Filter 4: bc found'] += 1
                    break
            else:
                # barcode not found
                stats['Filter 4: bc decode tried and failed'] += 1
                continue
    
            obs_umi, obs_tail_end = bc_umi_aligners[bc].get_umi_and_tail_end_pos(ops.seq_after_primer)
            umi_len_cntr[len(obs_umi)] += 1
            ops.set_barcode_and_umi(bc, obs_umi)
            ops.set_prefix_end_from_rel_prefix_end_after_primer(obs_tail_end)
            SeqIO.write(ops.rec_after_prefix, fq_out, 'fastq')
        stats['Total reads'] += i
        fq_out.close()

        for k, v in stats.items():
            cum_stats[k] += v
        for k, v in umi_len_cntr.items():
            cum_umi_len_cntr[k] += v
    
        if len(arguments.fastq_files) > 1:
            out_stats_fpath = base_out_fpath + '.stats.txt'
            write_stats_file_from_cntr(stats, out_stats_fpath)
            out_umi_fpath = base_out_fpath + '.umi_lens.txt'
            write_stats_file_from_cntr(umi_len_cntr, out_umi_fpath)

    out_stats_fpath = os.path.join(arguments.output_dir, 'cumulative.stats.txt')
    write_stats_file_from_cntr(cum_stats, out_stats_fpath)

    out_umi_fpath = os.path.join(arguments.output_dir, 'cumulative.umi_lens.txt')
    write_stats_file_from_cntr(cum_umi_len_cntr, out_umi_fpath)

    out_umi_fpath = os.path.join(arguments.output_dir, 'cumulative.umi_lens.pdf')
    fig, ax = umi_len_plot(cum_umi_len_cntr)
    fig.savefig(out_umi_fpath)
    

def demultiplex_bcs(fwd_primer_max_end, rc_primer_max_end, bc_oi_list, arguments):
    def make_base_out_fpath(fpath):
        fname = os.path.basename(fpath)
        for ending in ['.fastq', '.fq', '.fastq.gz', '.fq.gz']:
            if fname.endswith(ending):
                fname = fname[:-len(ending)] + '.oriented_and_demult.d{}+{}'.format(
                        arguments.max_err_decode,
                        arguments.reject_delta
                        )
                return os.path.join(arguments.output_dir, fname)
        raise ValueError(f'Unrecognized fastq ending in {fpath}')

    log.info('Loading barcode decoder...')
    bd = freebarcodes.decode.FreeDivBarcodeDecoder()
    bd.build_codebook_from_random_codewords(bc_oi_list, arguments.max_err_decode, arguments.reject_delta)

    log.info('Demultiplexing cell barcodes...')
    cum_stats = Counter()
    for fastq_fpath in arguments.fastq_files:
        log.info(f'Demultiplexing {fastq_fpath}')
        base_out_fpath = make_base_out_fpath(fastq_fpath)
        out_fastq_fpath = base_out_fpath + '.fastq'
        fq_out = open(out_fastq_fpath, 'w')
        log.info(f'Writing to {out_fastq_fpath}')
        stats = Counter()
        for i, rec in enumerate(SeqIO.parse(gzip_friendly_open(fastq_fpath), 'fastq')):
            if i % 1000000 == 0:
                log.debug(f'\t{i:,d}')

            ops = OrientedPrimerSeq(rec, fwd_primer_max_end, rc_primer_max_end)
            if ops.orientation_failed:
                stats['Filter 1: orientation failed'] += 1
                continue
            stats['Filter 1: orientation successful'] += 1
            if not ops.is_well_formed:
                stats['Filter 2: primer alignment fails check'] += 1
                continue
            stats['Filter 2: primer alignment passes check'] += 1
            if len(ops.seq_after_primer) < 18:
                stats['Filter 3: seq too short'] += 1
                continue
            stats['Filter 3: seq long enough'] += 1
            for obs_bc, start_pos in ops.get_kmers_at_primer_end_plusminus(16, 2):
                bc = bd.decode(obs_bc)
                if isinstance(bc, str):
                    ops.primer_end = start_pos
                    stats['Filter 4: bc found'] += 1
                    break
            else:
                # barcode not found
                stats['Filter 4: bc decode tried and failed'] += 1
                continue
    
            obs_umi = 'N'*10
            obs_tail_end = 16
            ops.set_barcode_and_umi(bc, obs_umi)
            ops.set_prefix_end_from_rel_prefix_end_after_primer(obs_tail_end)
            SeqIO.write(ops.rec_after_prefix, fq_out, 'fastq')
        stats['Total reads'] += i
        fq_out.close()

        for k, v in stats.items():
            cum_stats[k] += v
        if len(arguments.fastq_files) > 1:
            out_stats_fpath = base_out_fpath + '.stats.txt'
            write_stats_file_from_cntr(stats, out_stats_fpath)
    
    out_stats_fpath = os.path.join(arguments.output_dir, 'cumulative.stats.txt')
    write_stats_file_from_cntr(cum_stats, out_stats_fpath)
