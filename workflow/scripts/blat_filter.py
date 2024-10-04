import sys
import os
import pysam
import numpy as np
import pandas as pd
import time
import subprocess
from collections import OrderedDict
import argparse


class BlatAligner:
    """
    A class representing a blat aligner subprocess.

    Attributes:
        blat: path to boolean executable
        ref_2bit: path to .2bit formatted reference
        ooc: path to over-occurring N-mers file

    """

    def __init__(self, blat_executable, ref_2bit, ooc):
        self.blat = blat_executable
        self.ref_2bit = ref_2bit
        self.ooc = ooc

    def query(self, batch_query_filename, output_filename):
        completedProcess = subprocess.run(
            [self.blat, self.ref_2bit, batch_query_filename, "-ooc={0}".format(self.ooc), output_filename],
            capture_output=True, check=True)
        return completedProcess.returncode, completedProcess.stdout, completedProcess.stderr


class InsertSizeMetrics:
    """
    Class representing insert size metrics table
    """

    def __init__(self, insert_size_metrics_fn, alpha):
        _, _, insert_size_distr_df = self.parse_file(insert_size_metrics_fn)

        # get alpha/2 -> 1-alpha/2 percentile interval of insert sizes
        cumfrac = insert_size_distr_df['counts'].cumsum()/insert_size_distr_df['counts'].sum()
        self.lower_insert_size, self.upper_insert_size = cumfrac.loc[(cumfrac >= alpha/2) & (cumfrac <= (1 - alpha/2))].iloc[[0, -1]].index

    @staticmethod
    def parse_file(fn):

        with open(fn, 'r') as fh:
            line = fh.readline()
            while 'METRICS CLASS' not in line:
                line = fh.readline()

            metrics_cols_line = fh.readline()
            metrics_cols = metrics_cols_line.split('\t')
            metrics_values_line = fh.readline()
            metrics_values = metrics_values_line.split('\t')
            metrics_dict = {metrics_cols[i]: metrics_values[i] for i in range(len(metrics_cols))}

            width_percent_dict = {int(col.split('_')[2]): int(metrics_dict[col]) for col in
                                  [c for c in metrics_cols if 'WIDTH_OF_' in c]}

            while 'HISTOGRAM' not in line:
                line = fh.readline()

            histogram_cols_line = fh.readline()
            insert_sizes = []
            counts = []
            line = fh.readline()
            while line.strip():
                insert_size, count = line.split('\t')[:2] # CollectInsertSizeMetrics reports TANDEM and FR counts in additional columns if they occur in more than 5% of reads by default. Keep only first count column.
                insert_size = int(insert_size)
                count = int(count)

                insert_sizes.append(insert_size)
                counts.append(count)
                line = fh.readline()

        insert_size_distr_df = pd.DataFrame({'insert_size': insert_sizes, 'counts': counts})
        insert_size_distr_df = insert_size_distr_df.set_index('insert_size')

        return metrics_dict, width_percent_dict, insert_size_distr_df


class InputMAF:
    """
    A class representing contents of input MAF file.

    Attributes:
        input_maf_file: path to maf file
        variant_df: pandas dataframe listing all of variants in maf file along with their ref and alt counts
    """

    def __init__(
            self,
            maf,
            set_maf_col_manually=False,
            maf_chr_col=None,
            maf_start_pos_col=None,
            maf_end_pos_col=None,
            maf_ref_allele_col=None,
            maf_alt_allele_col=None,
            maf_t_ref_count_col=None,
            maf_t_alt_count_col=None,
            maf_chr_col_to_score_reads=None,
            maf_start_pos_col_to_score_reads=None,
            maf_pos_col_to_score_reads_empty_value=None
    ):
        """
        Parameters
        ----------
        maf_chr_col_to_score_reads: str
            Column name in maf referring to chromosome from a liftover task to a target reference for scoring reads.
        maf_start_pos_col_to_score_reads: str
            Column name in maf referring to start position from a liftover task to a target reference for scoring reads.
        maf_pos_col_to_score_reads_empty_value: Any
            Value indicating to skip the mutation because there is not a valid liftover coordinate in target reference
        """
        self.input_maf_file = maf
        input_maf = open(maf, "r")

        header = input_maf.readline()
        while header[0] == "#" or not header.strip():
            header = input_maf.readline()

        self.header = header.strip("\n").split("\t")

        if set_maf_col_manually:
            if None in [
                maf_chr_col,
                maf_start_pos_col,
                maf_end_pos_col,
                maf_ref_allele_col,
                maf_alt_allele_col,
                maf_t_ref_count_col,
                maf_t_alt_count_col,
                maf_chr_col_to_score_reads,
                maf_start_pos_col_to_score_reads,
                maf_pos_col_to_score_reads_empty_value
            ]:
                raise ValueError(
                    'set_maf_col_manually={0}, but one of the inputs in None.'.format(set_maf_col_manually)
                )
            self.contig_i = self.header.index(maf_chr_col)
            self.start_pos_i = self.header.index(maf_start_pos_col)
            self.end_pos_i = self.header.index(maf_end_pos_col)
            self.ref_allele_i = self.header.index(maf_ref_allele_col)
            self.alt_allele_i = self.header.index(maf_alt_allele_col)
            self._t_ref_count_i = self.header.index(maf_t_ref_count_col)
            self._t_alt_count_i = self.header.index(maf_t_alt_count_col)

            self.contig_to_score_reads_i = self.header.index(maf_chr_col_to_score_reads)
            self.start_pos_to_score_reads_i = self.header.index(maf_start_pos_col_to_score_reads)
        else:
            try:  # maflite
                self.contig_i = self.header.index("chr")
                self.start_pos_i = self.header.index("start")
                self.end_pos_i = self.header.index("end")
                self.ref_allele_i = self.header.index("ref_allele")
                self.alt_allele_i = self.header.index("tum_allele2")
            except:  # maf - this is what funcotate produces
                self.contig_i = self.header.index("Chromosome")
                self.start_pos_i = self.header.index("Start_position")
                self.end_pos_i = self.header.index("End_position")
                self.ref_allele_i = self.header.index("Reference_Allele")
                self.alt_allele_i = self.header.index("Tumor_Seq_Allele2")

            self.contig_to_score_reads_i = self.contig_i
            self.start_pos_to_score_reads_i = self.start_pos_i
            self.end_pos_to_score_reads_i = self.end_pos_i

            try:
                self._t_ref_count_i = self.header.index("t_ref_count")
                self._t_alt_count_i = self.header.index("t_alt_count")
            except:  # abs maf
                self._t_ref_count_i = self.header.index("ref_cnt")
                self._t_alt_count_i = self.header.index("alt_cnt")

        self.to_score_reads_empty_value = maf_pos_col_to_score_reads_empty_value

        self._maf_dict = {}
        variant_list = []

        for idx, line in enumerate(input_maf):
            maf_values = line.strip("\n").split("\t")

            try:
                ref_count = int(float(maf_values[self._t_ref_count_i]))
                alt_count = int(
                    float(maf_values[self._t_alt_count_i]))  # PCAWG mafs have alt/ref counts as 6.0 for some reason.
            except ValueError as ve:
                print("Detected mutation with null ref or alt count: {0}".format(line))
                raise ve

            key = '_'.join(
                [maf_values[0], maf_values[self.contig_i], maf_values[self.start_pos_i], maf_values[self.end_pos_i],
                 maf_values[self.ref_allele_i], maf_values[self.alt_allele_i]])
            if key in self._maf_dict:
                print("Detected multiple entries in MAF with same key:{0}\nline:{1}".format(key, line))
                raise ValueError("duplicate entries in MAF")

            self._maf_dict[key] = {'maf_record_values': maf_values, 'ref_count': ref_count, 'alt_count': alt_count}
            variant_list.append([key, ref_count, alt_count])

        df = pd.DataFrame(variant_list, columns=['variant', 'ref_count', 'alt_count'])
        df.set_index('variant', inplace=True)
        self.variant_df = df


class RealignmentFilter:
    """
    A class representing the blat realignment filter.

    """

    @staticmethod
    def score_read(chrm, pos, blat_alignments):

        """Score a read relative to intended genomic position."""

        # AS = blat alignment score, which for now is the number of matches in the alignment
        highest_AS_covering_mutation = 0
        highest_AS_not_covering_mutation = 0
        best_alignment = blat_alignments[0]
        for alignment in blat_alignments:

            alignment_contig = alignment["tName"]
            alignment_start = int(alignment["tStart"])
            alignment_end = int(alignment["tEnd"])
            alignment_matches = int(alignment["matches"])

            if ((alignment_contig == chrm) and (alignment_start < pos) and (alignment_end >= pos)):
                highest_AS_covering_mutation = max(alignment_matches, highest_AS_covering_mutation)

                if alignment_matches > highest_AS_covering_mutation:
                    best_alignment = alignment
            else:
                highest_AS_not_covering_mutation = max(alignment_matches, highest_AS_not_covering_mutation)

        primary_align_contig = best_alignment['tName']
        primary_align_start_pos = int(best_alignment["tStart"])
        primary_align_end_pos = int(best_alignment["tEnd"])
        return highest_AS_covering_mutation - highest_AS_not_covering_mutation, primary_align_contig, primary_align_start_pos, primary_align_end_pos

    @staticmethod
    def score_read_mate(primary_contig, primary_start_pos, primary_end_pos, primary_is_forward, blat_alignments,
                        lower_insert_size, upper_insert_size):
        highest_AS_covering_mutation = 0
        highest_AS_not_covering_mutation = 0
        for alignment in blat_alignments:

            alignment_contig = alignment["tName"]
            alignment_start = int(alignment["tStart"])
            alignment_end = int(alignment["tEnd"])
            alignment_matches = int(alignment["matches"])

            insert_size = (alignment_end - primary_start_pos) if primary_is_forward else (
                        primary_end_pos - alignment_start)
            valid_insert_size = (primary_contig == alignment_contig) and \
                                (insert_size > lower_insert_size) and \
                                (insert_size < upper_insert_size)

            if valid_insert_size:
                highest_AS_covering_mutation = max(alignment_matches, highest_AS_covering_mutation)
            else:
                highest_AS_not_covering_mutation = max(alignment_matches, highest_AS_not_covering_mutation)

        return highest_AS_covering_mutation - highest_AS_not_covering_mutation

    def __init__(
            self,
            blat_executable,
            bam,
            maf,
            database,
            ooc,
            paired_ended = True,
            insert_size_metrics = None,
            out_name=None,
            mapping_quality_threshold=5,
            base_quality_threshold=5,
            diff_score_concordant_threshold=3,
            diff_score_discordant_threshold=-10,
            concordant_fraction_threshold=0.75,
            set_maf_col_manually=False,
            maf_chr_col=None,
            maf_start_pos_col=None,
            maf_end_pos_col=None,
            maf_ref_allele_col=None,
            maf_alt_allele_col=None,
            maf_t_ref_count_col=None,
            maf_t_alt_count_col=None,
            maf_chr_col_to_score_reads=None,
            maf_start_pos_col_to_score_reads=None,
            maf_pos_col_to_score_reads_empty_value=None,
            mate_insert_size_alpha=0.05,
    ):

        """
        Parameters
        ----------
        blat_executable : str
            path to blat executable file
        bam : str
            path to tumor bam
        maf : str
            path to maf
        database : str
            path to .2bit formatted genome reference
        ooc : str
            path to over-ocurring N-mer file
        paired_ended: bool
            whether this is a paired end BAM. default True
        insert_size_metrics: str
            path to CollectMultipleMetrics insert size file from picard.analysis.InsertSizeMetrics
        out_name : str, optional
            base name for output mafs
        mapping_quality_threshold : int, optional
            mapping quality threshold (default=5) used to disqualify reads (ignore reads with mapping quality <= threshold)
        base_quality_threshold : int, optional
            base quality threshold (default=5) at variant site used to disqualify reads (ignores reads with base quality at variant position <= threshold)
        diff_score_concordant_threshold : int, optional
            scoring threshold (default=5) for classifying a read or mate as concordant (concordant if score > threshold)
        diff_score_discordant threshold : int optional
            scoring threshold (default=-10) for classifying a read as discordant (discordant if score < threshold)
        concordant_fraction_theshold : float, optional
            threshold (default = 0.75) for fraction of reads that are concordant for filtering of variants (PASS if concordant read fraction > threshold)

        # For liftover coordinates
        set_maf_col_manually: bool
            To set the maf columns manually or not. See class InputMAF for remaining arguments
                - maf_chr_col=None,
                - maf_start_pos_col=None,
                - maf_end_pos_col=None,
                - maf_ref_allele_col=None,
                - maf_alt_allele_col=None,
                - maf_t_ref_count_col=None,
                - maf_t_alt_count_col=None,
                - maf_chr_col_to_score_reads=None,
                - maf_start_pos_col_to_score_reads=None,
                - maf_pos_col_to_score_reads_empty_value=None

        mate_insert_size_alpha: float
            threshold for p-value test on insert size given the alignment of the mate
        """

        self._input_bam = pysam.AlignmentFile(bam, "rb")
        self._blat_aligner = BlatAligner(blat_executable, database, ooc)
        self._input_maf = InputMAF(
            maf,
            set_maf_col_manually=set_maf_col_manually,
            maf_chr_col=maf_chr_col,
            maf_start_pos_col=maf_start_pos_col,
            maf_end_pos_col=maf_end_pos_col,
            maf_ref_allele_col=maf_ref_allele_col,
            maf_alt_allele_col=maf_alt_allele_col,
            maf_t_ref_count_col=maf_t_ref_count_col,
            maf_t_alt_count_col=maf_t_alt_count_col,
            maf_chr_col_to_score_reads=maf_chr_col_to_score_reads,
            maf_start_pos_col_to_score_reads=maf_start_pos_col_to_score_reads,
            maf_pos_col_to_score_reads_empty_value=maf_pos_col_to_score_reads_empty_value
        )
        if not out_name:
            self._out_name = maf.rpartition("/")[-1].rpartition(".maf")[0]
        else:
            self._out_name = out_name

        self.paired_ended = paired_ended
        self._insert_size_metrics = InsertSizeMetrics(insert_size_metrics, float(mate_insert_size_alpha)) \
          if self.paired_ended else None
        self._mapping_quality_threshold = int(mapping_quality_threshold)
        self._base_quality_threshold = int(base_quality_threshold)
        self._diff_score_concordant_threshold = float(diff_score_concordant_threshold)
        self._diff_score_discordant_threshold = float(diff_score_discordant_threshold)
        self._concordant_fraction_threshold = float(concordant_fraction_threshold)
        self._blat_read_query_filename = out_name + '_query.fa'
        self._blat_read_query_output_psl_filename = out_name + '_output.psl'
        self._blat_mate_query_filename = out_name + '_mate_query.fa'
        self._blat_mate_query_output_psl_filename = out_name + '_mate_output.psl'


    def run_blat_aligner(self, mate=False):

        """Run blat aligner on a batch of queries."""

        if not mate:
            query_filename = self._blat_read_query_filename
            psl_output_filename = self._blat_read_query_output_psl_filename
        else:
            query_filename = self._blat_mate_query_filename
            psl_output_filename = self._blat_mate_query_output_psl_filename

        returncode, stdout, stderr = self._blat_aligner.query(query_filename, psl_output_filename)
        print('blat aligner return code = {}'.format(returncode))

    def collect_blat_queries(self):

        """Collect reads to be realigned in first pass."""

        query_list = []

        with open(self._blat_read_query_filename, 'w') as blat_query_file:

            # get pileups for each variant
            for key, value in self._input_maf._maf_dict.items():

                print("maf record key = {0}".format(key))

                contig_for_blat = "MT" if value['maf_record_values'][self._input_maf.contig_i] == "M" else \
                value['maf_record_values'][self._input_maf.contig_i]  ##fix inconsistency
                start_pos = int(value['maf_record_values'][self._input_maf.start_pos_i])

                contig_for_blat_to_score_reads = "MT" if value['maf_record_values'][
                                                             self._input_maf.contig_to_score_reads_i] == "M" else \
                value['maf_record_values'][self._input_maf.contig_to_score_reads_i]  ##fix inconsistency
                start_pos_to_score_reads = int(value['maf_record_values'][self._input_maf.start_pos_to_score_reads_i])

                if start_pos_to_score_reads == self._input_maf.to_score_reads_empty_value:
                    continue

                alt_allele = value['maf_record_values'][self._input_maf.alt_allele_i]
                ref_allele = value['maf_record_values'][self._input_maf.ref_allele_i]
                alt_count = value['alt_count']
                depth = value['alt_count'] + value['ref_count']

                read_names = set()

                # Chip note: do we check flags if a read is a supplemental alignment - feature in BWMem - dealing with split reads?
                #            Mutect1 does not consider split reads...so should not consider.

                for pileupcolumn in self._input_bam.pileup(contig_for_blat, start_pos - 1, start_pos, max_depth=2000000,
                                                           multiple_iterators=True):
                    if pileupcolumn.pos != start_pos - 1:
                        continue

                    # prcoess each read in the pileup
                    for read in pileupcolumn.pileups:

                        # READ FILTERING

                        # skip if read if satisfies any of a series of conditions
                        if read.alignment.query_name in read_names or read.alignment.is_duplicate or read.alignment.mapping_quality <= self._mapping_quality_threshold or read.alignment.is_qcfail or read.alignment.is_secondary or read.alignment.is_supplementary or read.alignment.is_unmapped:
                            continue

                        # set base_score_filter to True if  base quality score at query position is <=5
                        # note: read.query_position can be none for deletions
                        base_score_filter = (read.alignment.query_qualities[read.query_position] <= self._base_quality_threshold) if read.query_position is not None else False
                        # if read.query_position is None, but mutation is not a deletion, skip because it is not the alt we want.
                        base_score_filter = True if read.query_position is None and alt_allele != "-" else base_score_filter

                        # skip if base_score_filter == True
                        if base_score_filter:
                            continue

                        # DETERMINE WHETHER READ MAPPED as alt or indel at query base

                        is_alt = False
                        if alt_allele != "-" and ref_allele != "-":  # indels handled below
                            is_alt = True if read.alignment.seq[read.query_position] == alt_allele else False

                        elif len(read.alignment.blocks) > 1:
                            # is indel

                            blk_prev = read.alignment.blocks[0][1]
                            for blck in read.alignment.blocks[1:]:
                                blck_size = blck[0] - blk_prev
                                if blck_size == 0:  # insertion
                                    if blk_prev == start_pos:
                                        is_alt = True
                                else:  # deletion
                                    if blk_prev <= start_pos - 1 <= blk_prev + blck_size:
                                        is_alt = True

                                blk_prev = blck[1]

                        # add to set of query names
                        read_names.add(read.alignment.query_name)

                        if is_alt:

                            expanded_qname = read.alignment.query_name + ':' + key

                            query_list.append(
                                [expanded_qname, key, contig_for_blat_to_score_reads, start_pos_to_score_reads, None])

                            blat_query_file.write(">{0}\n".format(expanded_qname))
                            blat_query_file.write("{0}\n".format(read.alignment.query_sequence))

                        else:
                            pass

        df = pd.DataFrame(query_list, columns=['qname', 'variant', 'contig', 'pos', 'read_score']).astype(
            {'read_score': float})
        df.set_index('qname', inplace=True)
        self._supporting_reads_df = df

    # note: currently scoring is based on match counts returned by blat; we may want to consider basing scores on blat alignment scores
    # that would include gap penalties
    def score_supporting_reads(self, mate=False):

        """Score all of the realigned reads based on alignments returned by blat."""

        if not mate:
            blat_output_filename = self._blat_read_query_output_psl_filename
        else:
            blat_output_filename = self._blat_mate_query_output_psl_filename

        psl_fields = ['matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert',
                      'tBaseInsert',
                      'strand', 'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd',
                      'blockCount', 'blockSizes', 'qStarts', 'tStarts']
        matches_index = psl_fields.index('matches')
        qName_index = psl_fields.index('qName')
        tName_index = psl_fields.index('tName')
        tStart_index = psl_fields.index('tStart')
        tEnd_index = psl_fields.index('tEnd')

        with open(blat_output_filename, 'r') as psl_file:
            next_line = psl_file.readline()
            while not next_line.startswith('---'):
                next_line = psl_file.readline()

            last_query_name = None
            blat_alignment_list = []
            for psl_row in psl_file:
                psl_list = psl_row.split('\t')
                query_name = psl_list[qName_index]
                blat_alignment = {"matches": int(psl_list[matches_index]),
                                  "tName": psl_list[tName_index],
                                  "tStart": int(psl_list[tStart_index]),
                                  "tEnd": int(psl_list[tEnd_index])}
                if query_name == last_query_name:
                    blat_alignment_list.append(blat_alignment)
                else:
                    if last_query_name is not None:
                        if not mate:
                            query_row = self._supporting_reads_df.loc[last_query_name]
                            score, primary_align_contig, primary_align_start_pos, primary_align_end_pos = RealignmentFilter.score_read(
                                query_row['contig'], query_row['pos'], blat_alignment_list)
                            self._supporting_reads_df.loc[last_query_name, ['read_score']] = [score]
                            self._supporting_reads_df.loc[last_query_name, ['primary_align_contig']] = [
                                primary_align_contig]
                            self._supporting_reads_df.loc[last_query_name, ['primary_align_start_pos']] = [
                                primary_align_start_pos]
                            self._supporting_reads_df.loc[last_query_name, ['primary_align_end_pos']] = [
                                primary_align_end_pos]
                        else:
                            query_row = self._ambiguous_supporting_reads_df.loc[last_query_name]

                            # score = RealignmentFilter.score_read(query_row['contig'], query_row['mate_pos'], blat_alignment_list)
                            score = RealignmentFilter.score_read_mate(
                                primary_contig=query_row['primary_align_contig'].item(),
                                primary_start_pos=query_row['primary_align_start_pos'].item(),
                                primary_end_pos=query_row['primary_align_end_pos'].item(),
                                primary_is_forward=query_row['primary_align_is_forward'].item(),
                                blat_alignments=blat_alignment_list,
                                lower_insert_size=self._insert_size_metrics.lower_insert_size,
                                upper_insert_size=self._insert_size_metrics.upper_insert_size
                            )
                            self._ambiguous_supporting_reads_df.loc[last_query_name, ['mate_score']] = [score]
                    last_query_name = query_name
                    blat_alignment_list = [blat_alignment]

            # complete processing of last supporting read
            if last_query_name is not None:
                if not mate:
                    query_row = self._supporting_reads_df.loc[last_query_name]

                    score, primary_align_contig, primary_align_start_pos, primary_align_end_pos = RealignmentFilter.score_read(
                        query_row['contig'], query_row['pos'], blat_alignment_list)
                    self._supporting_reads_df.loc[last_query_name, ['read_score']] = [score]
                    self._supporting_reads_df.loc[last_query_name, ['primary_align_contig']] = [primary_align_contig]
                    self._supporting_reads_df.loc[last_query_name, ['primary_align_start_pos']] = [primary_align_start_pos]
                    self._supporting_reads_df.loc[last_query_name, ['primary_align_end_pos']] = [primary_align_end_pos]
                else:
                    query_row = self._ambiguous_supporting_reads_df.loc[last_query_name]
                    score = RealignmentFilter.score_read_mate(
                        primary_contig=query_row['primary_align_contig'].item(),
                        primary_start_pos=query_row['primary_align_start_pos'].item(),
                        primary_end_pos=query_row['primary_align_end_pos'].item(),
                        primary_is_forward=query_row['primary_align_is_forward'].item(),
                        blat_alignments=blat_alignment_list,
                        lower_insert_size=self._insert_size_metrics.lower_insert_size,
                        upper_insert_size=self._insert_size_metrics.upper_insert_size
                    )
                    self._ambiguous_supporting_reads_df.loc[last_query_name, ['mate_score']] = [score]

    def first_pass_process_read_scores(self):

        """Classify supporting reads (concordant, discordant, ambiguous) realigned in first pass."""

        # Here we classify reads as concordant, discordant or ambiguous
        # The classification algorithm should be reviewed and is easily modified.
        def classify_supporting_reads(row):

            read_score = row['read_score']
            concordant = 0
            discordant = 0
            ambiguous = 0

            if type(read_score) == float:
                if read_score > self._diff_score_concordant_threshold:
                    concordant = 1
                elif read_score < self._diff_score_discordant_threshold:
                    discordant = 1
                else:
                    ambiguous = 1
            else:
                discordant = 1

            return pd.Series([concordant, discordant, ambiguous], index=['concordant', 'discordant', 'ambiguous'])

        self._supporting_reads_df = self._supporting_reads_df.join(
            self._supporting_reads_df.apply(classify_supporting_reads, axis=1))

        self._variant_df = self._input_maf.variant_df.join(
            self._supporting_reads_df.groupby('variant')[['concordant', 'discordant', 'ambiguous']].sum())
        self._variant_df = self._variant_df.join(
            self._supporting_reads_df.groupby('variant').size().rename('read_count'))
        self._variant_df = self._variant_df.join(
            pd.DataFrame({'read_scores': self._supporting_reads_df.groupby('variant')['read_score'].apply(list)}))

    def first_pass_filter_judgements(self):

        """Judge variant (PASS, REJECT, or UNCERTAIN) based on classifications of supporting reads that were realigned in first pass."""

        def judge_variant(row):
            num_discordant = row['discordant']
            num_concordant = row['concordant']
            num_ambiguous = row['ambiguous']
            read_count = row['read_count']

            if num_concordant + num_ambiguous < 3:
                judgement = 'REJECT'
            elif num_concordant / read_count > self._concordant_fraction_threshold:
                judgement = 'PASS'
            elif num_discordant / read_count > 1 - self._concordant_fraction_threshold:
                judgement = 'REJECT'
            else:
                judgement = 'UNCERTAIN'

            return pd.Series([judgement], index=['judgement'])

        self._variant_df = self._variant_df.join(self._variant_df.apply(judge_variant, axis=1))

    def collect_mate_blat_queries(self):

        """Collect read mates to be realigned in second pass."""

        ambiguous_scoring_reads = self._supporting_reads_df.loc[self._supporting_reads_df['ambiguous'] == 1.0]
        ambiguous_scoring_read_set = set(ambiguous_scoring_reads.index.to_list())
        uncertain_variant_set = set(self._variant_df.loc[self._variant_df['judgement'] == 'UNCERTAIN'].index.to_list())

        query_mate_list = []
        blat_mate_query_count = 0

        with open(self._blat_mate_query_filename, 'w') as blat_query_file:

            # get pileups for each variant
            for key, value in self._input_maf._maf_dict.items():

                # we only consider variants that were classified UNCERTAIN in the first pass
                if key not in uncertain_variant_set:
                    continue

                print("maf record key = {0}".format(key))

                contig_for_blat = "MT" if value['maf_record_values'][self._input_maf.contig_i] == "M" else \
                value['maf_record_values'][self._input_maf.contig_i]
                start_pos = int(value['maf_record_values'][self._input_maf.start_pos_i])

                contig_for_blat_to_score_reads = "MT" if value['maf_record_values'][
                                                             self._input_maf.contig_to_score_reads_i] == "M" else \
                value['maf_record_values'][self._input_maf.contig_to_score_reads_i]  ##fix inconsistency
                start_pos_to_score_reads = int(value['maf_record_values'][self._input_maf.start_pos_to_score_reads_i])

                if start_pos_to_score_reads == self._input_maf.to_score_reads_empty_value:
                    continue

                read_names = set()

                # note that, by default, pileup will ignore orphaned reads (reads not in proper pair)
                for pileupcolumn in self._input_bam.pileup(contig_for_blat, start_pos - 1, start_pos, max_depth=2000000,
                                                           multiple_iterators=True):
                    if pileupcolumn.pos != start_pos - 1:
                        continue

                    # process each read in the pileup
                    for read in pileupcolumn.pileups:

                        # skip duplicates
                        if read.alignment.query_name in read_names:
                            continue
                        else:
                            read_names.add(read.alignment.query_name)

                        expanded_qname = read.alignment.query_name + ':' + key

                        # we only consider reads that were previously classified as ambiguous (neither clearly concordant nor discordant)
                        if expanded_qname not in ambiguous_scoring_read_set:
                            continue

                        primary_align_contig = ambiguous_scoring_reads.loc[expanded_qname, ['primary_align_contig']]
                        primary_align_start_pos = ambiguous_scoring_reads.loc[
                            expanded_qname, ['primary_align_start_pos']]
                        primary_align_end_pos = ambiguous_scoring_reads.loc[expanded_qname, ['primary_align_end_pos']]
                        primary_align_is_forward = not read.alignment.is_reverse

                        if read.alignment.is_proper_pair:
                            mate = self._input_bam.mate(read.alignment)

                            blat_query_file.write(">{0}\n".format(expanded_qname))
                            blat_query_file.write("{0}\n".format(mate.query_sequence))
                            blat_mate_query_count += 1
                            query_mate_list.append([
                                expanded_qname,
                                key,
                                contig_for_blat_to_score_reads,
                                start_pos_to_score_reads,
                                primary_align_contig,
                                primary_align_start_pos,
                                primary_align_end_pos,
                                primary_align_is_forward,
                                mate.reference_start + 2,
                                True,
                                None
                            ])

        df = pd.DataFrame(
            query_mate_list,
            columns=['qname', 'variant', 'contig', 'pos',
                     'primary_align_contig', 'primary_align_start_pos', 'primary_align_end_pos',
                     'primary_align_is_forward',
                     'mate_pos', 'mate_mapped', 'mate_score']
        ).astype({'mate_score': float})
        df.set_index('qname', inplace=True)
        self._ambiguous_supporting_reads_df = df
        return blat_mate_query_count

    def second_pass_process_read_scores(self):

        """Classify read mates (mate_concordant, mate_discordant) realigned in second pass."""

        def classify_ambiguous_supporting_reads(row):

            mate_score = row['mate_score']
            mate_concordant = 0
            mate_discordant = 0

            if type(mate_score) == float:
                if mate_score > self._diff_score_concordant_threshold:
                    mate_concordant = 1
                else:
                    mate_discordant = 1
            else:
                mate_discordant = 1

            return pd.Series([mate_concordant, mate_discordant], index=['mate_concordant', 'mate_discordant'])

        self._ambiguous_supporting_reads_df = self._ambiguous_supporting_reads_df.join(
            self._ambiguous_supporting_reads_df.apply(classify_ambiguous_supporting_reads, axis=1))

        self._uncertain_variant_df = self._variant_df.loc[self._variant_df['judgement'] == 'UNCERTAIN']
        self._uncertain_variant_df = self._uncertain_variant_df.join(
            self._ambiguous_supporting_reads_df.groupby('variant')[['mate_concordant', 'mate_discordant']].sum())
        self._uncertain_variant_df = self._uncertain_variant_df.join(
            self._ambiguous_supporting_reads_df.groupby('variant').size().rename('mate_count'))
        self._uncertain_variant_df = self._uncertain_variant_df.join(pd.DataFrame(
            {'mate_scores': self._ambiguous_supporting_reads_df.groupby('variant')['mate_score'].apply(list)}))

    def second_pass_filter_judgements(self):

        """Re-judge variants (PASS or REJECT) judged as UNCERTAIN in first pass, incorporating supporting mate classifications."""

        def judge_variant(row):
            num_discordant = row['discordant']
            num_concordant = row['concordant']
            num_mate_discordant = row['mate_discordant']
            num_mate_concordant = row['mate_concordant']
            read_count = row['read_count']

            # num_concordant + num_mate_concordant < 3, then reject
            if num_concordant + num_mate_concordant < 3:
                judgement = 'REJECT'
            elif (num_concordant + num_mate_concordant) / read_count > self._concordant_fraction_threshold:
                judgement = 'PASS'
            else:
                judgement = 'REJECT'

            return pd.Series([judgement], index=['judgement'])

        self._uncertain_variant_df = self._uncertain_variant_df.join(
            self._uncertain_variant_df.apply(judge_variant, axis=1), rsuffix='_second_pass')
        self._uncertain_variant_df['judgement'] = self._uncertain_variant_df['judgement_second_pass']
        self._uncertain_variant_df.drop('judgement_second_pass', axis=1, inplace=True)
        self._variant_df.update(self._uncertain_variant_df['judgement'])
        self._variant_df = self._variant_df.join(
            self._uncertain_variant_df[['mate_concordant', 'mate_discordant', 'mate_count', 'mate_scores']])

    def second_pass_reject_uncertain_variants(self):

        """Set judgement of UNCERTAIN variants from first pass to REJECT."""

        uncertain_variant_df = self._variant_df.loc[self._variant_df['judgement'] == 'UNCERTAIN']
        uncertain_variant_df['judgement'] = 'REJECT'
        self._variant_df.update(uncertain_variant_df['judgement'])
        self._variant_df['mate_count'] = np.nan
        self._variant_df['mate_concordant'] = np.nan
        self._variant_df['mate_discordant'] = np.nan

    def write_new_mafs(self):

        """Write three new mafs: new_maf contains both passed and rejected variants; debug_maf contains rejected variants; pass_maf contains passed variants."""

        with open(self._out_name + ".blat.all.maf", "w") as new_maf, open(self._out_name + ".blat.passed.maf",
                                                                          "w") as pass_maf, open(
                self._out_name + ".blat.rejected.maf", "w") as debug_maf:
            # write headers for three new mafs, including additional columns for blat filtering
            BLAT_FILTERING_COLUMNS = ["blat_read_concordant_count", "blat_read_discordant_count",
                                      "blat_read_ambiguous_count", "blat_read_scores",
                                      "blat_mate_concordant_count", "blat_mate_discordant_count", "blat_mate_scores",
                                      "FILTER_JUDGEMENT"]
            header = self._input_maf.header

            new_maf.write("\t".join(header + BLAT_FILTERING_COLUMNS) + "\n")
            pass_maf.write("\t".join(header + BLAT_FILTERING_COLUMNS) + "\n")
            debug_maf.write("\t".join(header + BLAT_FILTERING_COLUMNS) + "\n")

            def write_maf_record(row):

                # get entry in maf_dict
                maf_variant = self._input_maf._maf_dict[row.name]

                def float_to_int_to_str(x):
                    x = 'NaN' if np.isnan(x) else str(int(x))
                    return x

                # note: we write into the MAF up to 50 read and mate scores

                if np.isnan(row['concordant']):
                    read_columns = ['', '', '', '']
                else:
                    read_columns = [str(row['concordant']), str(row['discordant']), str(row['ambiguous']),
                                    ';'.join(map(float_to_int_to_str, row['read_scores'][:50]))]

                if not self.paired_ended or np.isnan(row['mate_concordant']):
                    mate_columns = ['', '', '']
                else:
                    mate_columns = [str(row['mate_concordant']), str(row['mate_discordant']),
                                    ';'.join(map(float_to_int_to_str, row['mate_scores'][:50]))]

                new_row = '\t'.join(
                    maf_variant['maf_record_values'] + read_columns + mate_columns + [row['judgement']]) + '\n'

                new_maf.write(new_row)
                if row['judgement'] == 'REJECT':
                    debug_maf.write(new_row)
                else:
                    pass_maf.write(new_row)

            self._variant_df.apply(write_maf_record, axis=1)


if __name__ == "__main__":

    # parser = argparse.ArgumentParser(description="Filter mutations with Blat.")
    # parser.add_argument('--blat_executable', help='path to blat executable (jar file)', required=True)
    # parser.add_argument('--bam', help='bam file', required=True)
    # parser.add_argument('--maf', help='maf file. Lists mutations to filter.', required=True)
    # parser.add_argument('--database', help='Reference 2bit file to blat sequences to.', required=True)
    # parser.add_argument('--ooc', help='ooc specific to reference. Run blat -makeOoc', required=True)
    # parser.add_argument(
    #     '--insert_size_metrics',
    #     help='Output metrics file from CollectMultipleMetrics or picard CollectInsertSizeMetrics',
    #     required=False
    # )
    # parser.add_argument(
    #     '--paired_ended',
    #     action="store_true"
    # )
    # parser.add_argument('--out_name', help='Output prefix for output maf files', default=None)
    # parser.add_argument('--mapping_quality_threshold', type=int, default=5)
    # parser.add_argument('--base_quality_threshold', type=int, default=5)
    # parser.add_argument('--diff_score_concordant_threshold', type=int, default=3)
    # parser.add_argument('--diff_score_discordant_threshold', type=int, default=-10)
    # parser.add_argument('--concordant_fraction_threshold', type=float, default=0.75)
    # parser.add_argument('--set_maf_col_manually', action='store_true')
    # parser.set_defaults(set_maf_col_manually=False)

    # parser.add_argument('--maf_chr_col', default='Chromosome')
    # parser.add_argument('--maf_start_pos_col', default='Start_position')
    # parser.add_argument('--maf_end_pos_col', default='End_position')
    # parser.add_argument('--maf_ref_allele_col', default='Reference_Allele')
    # parser.add_argument('--maf_alt_allele_col', default='Tumor_Seq_Allele2')
    # parser.add_argument('--maf_t_ref_count_col', default='t_ref_count')
    # parser.add_argument('--maf_t_alt_count_col', default='t_alt_count')
    # parser.add_argument('--maf_chr_col_to_score_reads', default=None)
    # parser.add_argument('--maf_start_pos_col_to_score_reads', default=None)
    # parser.add_argument('--maf_pos_col_to_score_reads_empty_value', default=-1)
    # parser.add_argument('--mate_insert_size_alpha', default=0.05)
    # args = parser.parse_args()
    # args_dict = vars(args)
    # print(args)
    args_dict = {
        "blat_executable": snakemake.params.blat_binary,
        "bam": snakemake.input.bam,
        "maf": snakemake.input.maf,
        "database": snakemake.input.database,
        "ooc": snakemake.input.occ,
        "paired_ended": True,
        "out_name": snakemake.params.output_prefix,
        "insert_size_metrics": snakemake.input.insert_size,
        "mapping_quality_threshold": 5,
        "base_quality_threshold": 5,
        "diff_score_concordant_threshold":3,
        "diff_score_discordant_threshold":10,
        "concordant_fraction_threshold": 0.75,
        "set_maf_col_manually": True,
        "maf_chr_col": "Chromosome",
        "maf_ref_allele_col": 'Reference_Allele',
        "maf_alt_allele_col": "Tumor_Seq_Allele2",
        "maf_chr_col_to_score_reads": "Chromosome",
        "maf_start_pos_col": "Start_Position",
        "maf_start_pos_col_to_score_reads": "Start_Position",
        "maf_end_pos_col": "End_Position",
        "maf_t_ref_count_col": "t_ref_count",
        "maf_t_alt_count_col": "t_alt_count",
        "maf_pos_col_to_score_reads_empty_value": -1,
        "mate_insert_size_alpha": 0.05
    }
    print(args_dict)
    os.makedirs(os.path.split(args_dict["out_name"])[0], exist_ok=True)
    RF = RealignmentFilter(**args_dict)

    # FIRST PASS

    print("\n**** first pass collecting blat queries ****\n")
    RF.collect_blat_queries()

    print("\n**** first pass run blat aligner ****\n")
    RF.run_blat_aligner()

    print("\n**** first pass score supporting reads ****\n")
    RF.score_supporting_reads()

    print("\n**** first pass classify reads ****\n")
    RF.first_pass_process_read_scores()

    print("\n**** first pass filter variants ****\n")
    RF.first_pass_filter_judgements()

    # SECOND PASS
    if RF.paired_ended:
        print("\n**** second pass collecting blat mate queries ****\n")
        blat_mate_query_count = RF.collect_mate_blat_queries()

        if blat_mate_query_count > 0:

            print("\n**** second pass run blat aligner on mates ****\n")
            RF.run_blat_aligner(mate=True)

            print("\n**** second pass score mates ****\n")
            RF.score_supporting_reads(mate=True)

            print("\n**** second pass classify mates ****\n")
            RF.second_pass_process_read_scores()

            print("\n**** second pass filter variants ****\n")
            RF.second_pass_filter_judgements()

        else:
            print("\n**** second pass: no mate queries ****")
            print("**** second pass reject UNCERTAIN variants ****\n")
            RF.second_pass_reject_uncertain_variants()

    # WRITE RESULTS TO MAF FILES

    print("\n**** write new mafs ****\n")
    RF.write_new_mafs()
