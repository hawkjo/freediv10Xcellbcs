from Bio import Align
from Bio import SeqIO
from .misc import gzip_friendly_open
from .constants import sp1_10X_primer


class PrimerAligner:
    primerseq = sp1_10X_primer
    err_thresh = -int(0.25*len(primerseq))  # no more than 25% errors
    aligner = Align.PairwiseAligner()
    aligner.wildcard = 'N'
    aligner.match = 0
    aligner.mismatch = -1
    aligner.gap_score = -1
    aligner.target_left_gap_score = 0
    aligner.target_right_gap_score = 0

    def __init__(self, fwd_primer_max_end=None, rc_primer_max_end=None):
        self.fwd_primer_max_end = fwd_primer_max_end
        self.rc_primer_max_end = rc_primer_max_end

    def align_rec(self, rec):
        return self.aligner.align(self.primer_seq, str(rec.seq))

    def align_seq(self, seq):
        return self.aligner.align(self.primer_seq, seq)

    def iterate_recs_with_primer_pos(self, fastq_files):
        """
        Iterates a list of fastq files and yields rec, strand, primer_start, primer_end
        """
        for fastq_file in fastq_files:
            for rec in SeqIO.parse(gzip_friendly_open(fastq_file), 'fastq'):
                strand, start, end = self.get_primer_pos_in_rec(rec)
                yield rec, strand, start, end

    def get_primer_pos_in_rec(self, rec):
        """
        Finds the strand, start, and end of primer in given rec. Returns all None if not found.
        """
        best_fwd_alignment = self.align_rec(rec[:self.fwd_primer_max_end])[0]
        best_rc_alignment = self.align_rec(rec.reverse_complement()[:self.rc_primer_max_end])[0]
        if best_fwd_alignment.score >= self.err_thresh and best_rc_alignment.score < self.err_thresh:
            start = best_fwd_alignment.aligned[1][0][0]  
            end = best_fwd_alignment.aligned[1][-1][-1] 
            return '+', start, end
        elif best_fwd_alignment.score < self.err_thresh and best_rc_alignment.score >= self.err_thresh:
            start = best_rc_alignment.aligned[1][0][0]
            end = best_rc_alignment.aligned[1][-1][-1]
            return '-', start, end
        else:
            return None, None, None


class OrientedPrimerSeq:
    def __init__(self, rec, fwd_primer_max_end, rc_primer_max_end):
        self.primer_aligner = PrimerAligner(fwd_primer_max_end, rc_primer_max_end)
        self.strand, self.primer_start, self.primer_end = self.primer_aligner.get_primer_pos_in_rec(rec)
        if self.strand is None:
            self.orientation_failed = True
            return
        else:
            self.orientation_failed = False
            if self.strand == '-':
                rec = rec.reverse_complement()
                name = f'{rec.name}:rev_comp'
                self.rec.id = name
                self.rec.name = name
                self.rec.description = rec.description
            
        self.seq = str(self.rec.seq)
        self.bc = None
        self.umi = None

        max_primer_end = fwd_primer_max_end if self.strand == '+' else rc_primer_max_end
        self.is_well_formed = bool(self.primer_start > 0 and self.primer_end < max_primer_end)  # has leading and trailing gaps
        
    @property
    def observed_primer(self):
        return self.seq[self.primer_start:self.primer_end]
    
    @property
    def observed_prefix(self):
        return self.seq[self.primer_start:self.prefix_end]

    @property
    def seq_after_primer(self):
        return self.seq[self.primer_end:]
    
    @property
    def rec_after_primer(self):
        return self.rec[self.primer_end:]
    
    @property
    def seq_after_prefix(self):
        return self.seq[self.prefix_end:]
    
    @property
    def rec_after_prefix(self):
        return self.rec[self.prefix_end:]
    
    def set_barcode_and_umi(self, bc, umi):
        assert self.bc is None and self.umi is None, (self.bc, self.umi)
        self.bc = bc
        self.umi = umi
        name = f'{bc}_{umi}#{self.rec.id}'
        self.rec.id = name
        self.rec.name = name
        
    def set_prefix_end_from_rel_prefix_end_after_primer(self, rel_prefix_end_after_primer):
        self.prefix_end = self.primer_end + rel_prefix_end_after_primer
        
    def get_kmers_at_primer_end_plusminus(self, k, plusminus):
        """
        Gets kmers starting at primer_end plus-minus plusminus bp.
        """
        # Iterate out from the center: 0, +1, -1, +2, -2, etc.
        starts = [self.primer_end]
        for abs_plusminus in range(1, plusminus + 1):
            starts.extend([self.primer_end + abs_plusminus, self.primer_end - abs_plusminus])
        for start in starts:
            subseq = self.seq[start:start+k]
            if len(subseq) == k:
                yield subseq, start


class BcUmiTailAligner:
    tail_5p_kit_tso = 'TTTCTTATATGGG'
    tail_3p_kit_polyT = 'T'*len(tail_5p_kit_tso)
    
    aligner = Align.PairwiseAligner()
    aligner.wildcard = 'N'
    aligner.mismatch = -1
    aligner.gap_score = -1.1
    aligner.target_left_gap_score = -1.9
    aligner.query_left_gap_score = -1.9
    aligner.target_right_gap_score = 0
        
    def __init__(self, bc, umi_len, kit):
        if kit == '3p':
            self.tail = self.tail_3p_kit_polyT
        elif kit == '5p':
            self.tail = self.tail_5p_kit_tso
        else:
            raise ValueError('kit must be either 3p or 5p')
            
        self.bc = bc
        self.umi_len = umi_len
        self.prefixes = [self.bc, 'N'*self.umi_len, self.tail]
        self.full_prefix = ''.join(self.prefixes)
        self.bc_end = len(self.bc)
        self.tail_start = len(self.bc) + self.umi_len
        self.tail_end = len(self.full_prefix)
        self.query_points_of_interest = [len(self.bc), len(self.bc) + umi_len, len(self.full_prefix)]
        self.max_query_len = int(1.5*len(self.full_prefix))
        
    def find_key_boundaries(self, seq):
        alignment = self.aligner.align(self.full_prefix, seq[:self.max_query_len])[0]
        obs_bc_end, obs_tail_start, obs_tail_end = None, None, None
        for i in range(len(alignment.aligned[0])):
            tstart, tend = alignment.aligned[0][i]
            qstart, qend = alignment.aligned[1][i]
            
            if obs_bc_end is None:
                if tstart <= self.bc_end <= tend:
                    obs_bc_end = qstart + self.bc_end - tstart
                elif tstart > self.bc_end:
                    if i == 0:  # bizarre alignment. discard
                        return None
                    obs_bc_end = alignment.aligned[1][i-1][1]  # prev_qend
                    
            if obs_tail_start is None:
                if tstart <= self.tail_start < tend:
                    obs_tail_start = qstart + self.tail_start - tstart
                elif tstart > self.tail_start:
                    obs_tail_start = tstart
                else:
                    continue
                break
                
        tstart, tend = alignment.aligned[0][-1]
        qstart, qend = alignment.aligned[1][-1]
        obs_tail_end = qend + self.tail_end - tend
        
        return obs_bc_end, obs_tail_start, obs_tail_end
        
    def get_umi_and_tail_end_pos(self, seq):
        obs_bc_end, obs_tail_start, obs_tail_end = self.find_key_boundaries(seq)
        umi = seq[obs_bc_end:obs_tail_start]
        return umi, obs_tail_end
    
    def get_umi(self, seq):
        return self.get_umi_and_tail_end_pos(seq)[0]
    
    def get_10bp_umi(self, seq):
        obs_bc_end, obs_tail_start, obs_tail_end = self.find_key_boundaries(seq)
        umi = seq[obs_bc_end:obs_bc_end+10]
        return umi
    
    def get_umi_and_tail_start_pos(self, seq):
        obs_bc_end, obs_tail_start, obs_tail_end = self.find_key_boundaries(seq)
        umi = seq[obs_bc_end:obs_tail_start]
        return umi, obs_tail_start
    
    def get_observed_bc_umi_tail(self, seq):
        obs_bc_end, obs_tail_start, obs_tail_end = self.find_key_boundaries(seq)
        return seq[0:obs_bc_end], seq[obs_bc_end:obs_tail_start], seq[obs_tail_start:obs_tail_end]
