import os
import logging
import numpy as np
import freebarcodes
import freebarcodes.decode
from Bio import SeqIO
from collections import Counter
from .discovery import find_primer_dist_ends, discover_bcs
from .misc import load_bc_list, gzip_friendly_open, write_stats_file_from_cntr
from .aligners import OrientedPrimerSeq, BcUmiTailAligner
from .plotting import umi_len_plot


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
