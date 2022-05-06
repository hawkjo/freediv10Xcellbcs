from Bio import Align
from .constants import sp1_10X_primer


class PrimerAligner:
    primerseq = sp1_10X_primer
    aligner = Align.PairwiseAligner()
    aligner.wildcard = 'N'
    aligner.match = 0
    aligner.mismatch = -1
    aligner.gap_score = -1
    aligner.target_left_gap_score = 0
    aligner.target_right_gap_score = 0

    def __init__(self):
        pass

    def align_rec(self, rec):
        return self.aligner.align(self.primer_seq, str(rec.seq))

    def align_seq(self, seq)
        return self.aligner.align(self.primer_seq, seq)


class OrientedPrimerSeq:
    primerseq = 'CTACACGACGCTCTTCCGATCT'
    
    aligner = Align.PairwiseAligner()
    aligner.wildcard = 'N'
    aligner.mismatch = -1
    aligner.gap_score = -2.1
    aligner.target_left_gap_score = 0
    aligner.target_right_gap_score = 0

    def __init__(self, fq_rec, hit):
        self.hit = hit
        if self.hit.strand == '+':
            self.fq_rec = fq_rec
        else:
            self.fq_rec = fq_rec.reverse_complement()
            name = f'rc-{fq_rec.name}'
            self.fq_rec.id = name
            self.fq_rec.name = name
            self.fq_rec.description = fq_rec.description
            
        self.seq = str(self.fq_rec.seq)
        self.bc = None
        self.umi = None
        
        start_buff = self.hit.pstart + 5
        end_buff = len(self.primerseq) - self.hit.pend + 5
        self.subseq_start = max(0, self.hit.qstart - start_buff)
        self.subseq_end = self.hit.qend + end_buff
        self.subseq = self.seq[self.subseq_start:self.subseq_end]
        
        self.alignment = self.aligner.align(self.primerseq, self.subseq)[0]
        self.is_well_formed = bool(self.alignment.path[1][0] == 0 and self.alignment.path[-2][0] == len(self.primerseq))  # has leading and trailing gaps
        self.primer_start = self.subseq_start + self.alignment.path[1][1]
        self.primer_end = self.subseq_start + self.alignment.path[-2][1]  # Query position where trailing gaps begin
        
    @property
    def observed_primer(self):
        return self.seq[self.primer_start:self.primer_end]
    
    @property
    def observed_prefix(self):
        return self.seq[self.primer_start:self.tail_end]
    
    @property
    def seq_after_primer(self):
        return self.seq[self.primer_end:]
    
    @property
    def rec_after_primer(self):
        return self.fq_rec[self.primer_end:]
    
    @property
    def seq_after_prefix(self):
        return self.seq[self.prefix_end:]
    
    @property
    def rec_after_prefix(self):
        return self.fq_rec[self.prefix_end:]
    
    def set_barcode_and_umi(self, bc, umi):
        assert self.bc is None and self.umi is None, (self.bc, self.umi)
        self.bc = bc
        self.umi = umi
        name = f'{bc}_{umi}#{self.fq_rec.id}'
        self.fq_rec.id = name
        self.fq_rec.name = name
        
    def set_prefix_end_from_rel_tail_end_after_primer(self, tail_end_after_primer):
        self.prefix_end = self.primer_end + tail_end_after_primer
        
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
    tail_3p_kit_polyA = 'T'*15
    
    aligner = Align.PairwiseAligner()
    aligner.wildcard = 'N'
    aligner.mismatch = -1
    aligner.gap_score = -1.1
    aligner.target_left_gap_score = -1.9
    aligner.query_left_gap_score = -1.9
    aligner.target_right_gap_score = 0
        
    def __init__(self, bc, umi_len, kit):
        if kit == '3p':
            self.tail = self.tail_3p_kit_polyA
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
                
        assert obs_bc_end is not None, str(alignment)
        assert obs_tail_start is not None, str(alignment)
        
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
    
    def get_bc_umi_tail(self, seq):
        obs_bc_end, obs_tail_start, obs_tail_end = self.find_key_boundaries(seq)
        return seq[0:obs_bc_end], seq[obs_bc_end:obs_tail_start], seq[obs_tail_start:obs_tail_end]



