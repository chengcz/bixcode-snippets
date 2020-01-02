#!/usr/bin/env python
# -*- coding:utf-8 -*-

from .log import logger
from .file import Files
from .interval import AtomicInterval, Interval, intervals
from .sequence import Sequence
from collections import OrderedDict


class StructureInforError(Exception):
    def __init__(self, info):
        self.message = info


class Transcript(object):
    '''
    Transcript structure annotation information

    parameter
    ---------
    Tid:        sring, Transcript IDs
    chro:       sting, chromesome ID
    interval:   Interval, 0-based location
    strand:     string, choices from {-, +}
    exon:       Interval,
    cds:        Interval,
    utr:        Interval,
    infor:      dict, attr of gtf file
    '''
    __slots__ = (
        '_id', '_chrom', '_interval', '_start', '_end', '_strand',
        '_exon', '_cds', '_utr', '_attri', '_utr5', '_utr3', '_msg'
    )
    def __init__(self, Tid, chro, interval, strand, exon=None,
        cds=None, utr=None, infor=None, strict=True
    ):
        self._id = Tid
        self._chrom = chro
        self._interval = interval
        self._start = interval.lower
        self._end = interval.upper
        self._strand = strand
        self._exon = exon if (exon is not None) else Interval()
        self._cds = cds if (cds is not None) else Interval()
        self._utr = utr if (utr is not None) else Interval()
        self._attri = infor if infor else {}
        ### adjustment and validation
        self.__Interval_self_adjustment()
        self._utr5 = Interval(*[x for x in self._utr if x < self.cds])
        self._utr3 = Interval(*[x for x in self._utr if x > self.cds])
        self._msg = '\n'.join([
            '-'*50,
            '    Transcript ID: {}'.format(Tid),
            '    Location:      {}:{}-{}:{}'.format(
                chro, interval.lower, interval.upper, strand
            ),
            '    Exon:          {}'.format(self._exon),
            '    Cds:           {}'.format(self._cds),
            '    Utr:           {}'.format(self._utr),
            '    Attri:         {}'.format(self._attri.items()),
            '-'*50,
        ])
        try:
            self.__data_validation(strict)
        except StructureInforError as e:
            if strict:
                logger.error('\n{}'.format(self._msg))
                raise StructureInforError(e)
            logger.warning('\n{}'.format(self._msg))
        if self._strand == '-':
            self._utr5, self._utr3 = self._utr3, self._utr5

    def __str__(self):
        return '<Transcript object>\n{}'.format(self._msg)

    def __repr__(self):
        return self.__str__()

    def __Interval_self_adjustment(self):
        flag_exon = not self._exon.is_empty()
        flag_cds = not self._cds.is_empty()
        flag_utr = not self._utr.is_empty()
        if flag_exon and flag_cds:
            if not flag_utr:
                self._utr = self._exon - self._cds
        elif flag_exon and not flag_cds:
            if flag_utr:
                self._cds = self._exon - self._utr
            else:
                self._cds = Interval()
                self._utr = Interval()
        elif not flag_exon and flag_cds:
            if flag_utr:
                self._exon = self._cds | self._utr
            else:
                self._exon = self._cds
                self._utr = Interval()
        else:
            raise StructureInforError('Exon & Cds interval both is missing.')

    def __data_validation(self, strict=True):
        if self._strand not in ('-', '+'):
            raise StructureInforError('unsupported strand Symbol, select from (-, +)')
        if self._interval.is_empty():
            raise StructureInforError('Transcript interval is wrong.')
        if self._exon.is_empty():
            raise StructureInforError('Exon interval is wrong.')

        if self._exon.to_atomic() != self._interval:
            raise StructureInforError('Unmatched interval between Transcript with Exon.')
        if (not self._cds.is_empty()) or (not self._utr.is_empty()):
            if self._exon != (self._cds | self._utr):
                raise StructureInforError('Unmatched interval between Exon with CDS & UTR.')
        if (not self._cds.is_empty()) and (not self._utr.is_empty()):
            if not (self._cds & self._utr).is_empty():
                raise StructureInforError('CDS interval Overlap with UTR interval.')

    def __len__(self):
        return self.length

    @property
    def id(self):
        return self._id

    @property
    def chrom(self):
        return self._chrom

    @property
    def start(self):
        '''
        1-based
        '''
        return self._start + 1

    @property
    def end(self):
        return self._end

    @property
    def strand(self):
        return self._strand

    @property
    def exon(self):
        return self._exon

    @property
    def exon_count(self):
        return len(self._exon)

    @property
    def length(self):
        return self._exon.length

    @property
    def intron(self):
        TransRegion= AtomicInterval(self.exon.lower, self.exon.upper)
        # intron = self._interval - self._exon
        # intron = Interval(
        #     *[x for x in intron if x.overlaps(TransRegion, False)]
        # )
        intron = TransRegion - self.exon
        return intron

    @property
    def cds(self):
        return self._cds

    @property
    def cds_count(self):
        return len(self._cds)

    @property
    def length_cds(self):
        return self._cds.length

    @property
    def utr(self):
        return self._utr

    @property
    def utr5(self):
        return self._utr5

    @property
    def utr3(self):
        return self._utr3

    @property
    def name(self):
        return self._attri.get('transcript_name', self._id)

    @property
    def biotype(self):
        tmp = self._attri
        return tmp.get('transcript_type', tmp.get('transcript_biotype', tmp.get('biotype', "")))

    @property
    def gene_id(self):
        return self._attri.get('gene_id', self._id)

    @property
    def gene_name(self):
        return self._attri.get('gene_name', "")

    @property
    def gene_biotype(self):
        tmp = self._attri
        return tmp.get('gene_type', tmp.get('gene_biotype', tmp.get('biotype', "")))

    @property
    def bed(self):
        '''
        basic bed:
            chro, start, end, name, score, strand
        '''
        return (self._chrom, self._start, self._end, self._id, 0, self._strand)

    def del_version(self, sep='.'):
        '''
        remove version infor of transcript id, defaule version infor at '.' after
        '''
        tid = self._id
        if sep in tid:
            self._id = tid[:tid.rindex(sep)]
        gid = self.gene_id
        if sep in gid:
            self._attri['gene_id'] = gid[:gid.rindex(sep)]

    def id_modifer(self, func):
        '''
        modify transcript_id, gene_id

        parameter
        ---------
        func:     function, recommend use lambda
        '''
        self._id = func(self._id)
        self._attri['gene_id'] = func(self.gene_id)

    def to_gtf(self, fp):
        '''
        parameter
        ---------
        fp:    file handle for output standrand GTF file
        '''
        attr = self.__attri_of_gtfline()

        transcript = (
            self._chrom, '.', 'transcript', self._start + 1, self._end,
            '.', self._strand, '.', attr,
        )
        Record = []
        tmp = zip(
            ('exon', 'CDS', '5UTR', '3UTR'),
            (self._exon, self._cds, self._utr5, self._utr3)
        )
        for label, regions in tmp:
            for region in regions:
                if region.is_empty():
                    continue
                start, end = region.lower, region.upper
                line = (self._chrom, '.', label, start + 1, end, '.', self._strand, '.', attr)
                Record.append(line)
        order = False if self._strand == '+' else True
        Record = sorted(Record, key=lambda x: x[3], reverse=order)

        Record.insert(0, transcript)
        for feature in Record:
            fp.write(self.__list2str(feature))

    def __attri_of_gtfline(self):
        skipkeys = {
            'gene_id', 'gene_name', 'gene_type', 'gene_biotype',
            'transcript_id', 'transcript_name', 'transcript_type',
            'transcript_biotype', 'Parent', 'ID', 'gene', 'gbkey', 'Name',
        }
        keepkeys = set(self._attri.keys()) - skipkeys
        attr = ['transcript_id "{}"; gene_id "{}"; '.format(self.id, self.gene_id), ]
        keys = ['transcript_name', 'gene_name', 'transcript_type', 'gene_type']
        for key in keys:
            value = self._attri.get(key, None)
            if value:
                attr.append('{} "{}"; '.format(key, value))
        attr += ['{} "{}"; '.format(i, self._attri[i]) for i in keepkeys]
        return ''.join(attr)

    @staticmethod
    def __list2str(lst):
        return '\t'.join([str(i) for i in lst]) + '\n'

    def to_bed(self, fp):
        '''
        parameter
        ---------
        fp:    file handle for output 12 columns bed file
        '''
        if not self._cds.is_empty():
            cstart, cend = self._cds.lower, self._cds.upper
        else:
            cstart, cend = (self.end, self.end)
        exon_num = self.exon_count
        exon_len = ''.join(['{},'.format(len(x)) for x in self._exon])
        exon_start = ''.join(['{},'.format(x.lower - self._start) for x in self._exon])

        Record = (
            self._chrom, self._start, self._end, self._id, 0, self._strand,
            cstart, cend, 0, exon_num, exon_len, exon_start,
        )
        fp.write(self.__list2str(Record))

    def to_genePred(self, fp, refFlat=True):
        '''
        parameter
        ---------
        fp:     file handle for output GenePred file,
        refFlat:   bool
            whether create refFlat style genePred, {True, False}

        reference
        ---------
        url: https://genome.ucsc.edu/FAQ/FAQformat#format9
        '''
        if not self._cds.is_empty():
            cstart, cend = self._cds.lower, self._cds.upper
        else:
            cstart, cend = (self.end, self.end)
        exon_num = self.exon_count
        estart = ''.join(['{},'.format(i.lower) for i in self._exon])
        eend = ''.join(['{},'.format(i.upper) for i in self._exon])

        Record = [
            self.name, self._chrom, self._strand, self._start, self._end,
            cstart, cend, exon_num, estart, eend,
        ]
        if refFlat:
            Record.insert(0, self.gene_id)
        fp.write(self.__list2str(Record))

    @property
    def summary(self):
        '''
        return
        ------
        id, chro, start, end, strand, name, type, length, gene_id, gene_name, gene_type
        '''
        summary = (
            self.id, self.chrom, self.start, self.end, self.strand,
            self.name, self.biotype, self.length,
            self.gene_id, self.gene_name, self.gene_biotype
        )
        return summary

    def overlap_with(self, other, flag_strand=True):
        '''
        '''
        flag = False
        if isinstance(other, Transcript):
            if self.chrom != other.chro:
                return False
            if flag_strand and self.strand != other.strand:
                return False
            if (self._exon & other.exon).is_empty():
                return False
            return True
        else:
            return NotImplemented

    @staticmethod
    def __extract_seq(faidx, chro, start, end):
        try:
            seq = faidx[chro][start:end]
        except KeyError:
            seq = ''
            logger.error('Chromesome("{}") does not exist in fasta file, Skip...'.format(chro))
        return seq

    def extract_genomic_seq(self, seqfp):
        '''
        parameter
        ---------
        seqfp: dict
            {chro: seq, }

        return Sequence object
        '''
        seq = self.__extract_seq(seqfp, self._chrom, self._start, self._end)
        obj = Sequence(
            '{}'.format(self._id),
            seq,
            'genomic sequence, gene_id:{}'.format(self.gene_id)
        )
        if self._strand == '-':
            obj = obj.reverse_complement()
        return obj

    def extract_transcript_seq(self, seqfp):
        '''
        return Sequence object
        '''
        seq = ''.join([
            self.__extract_seq(seqfp, self._chrom, x.lower, x.upper) for x in self._exon
        ])
        obj = Sequence(
            self._id, seq, 'transcript sequence, gene_id:{}'.format(self.gene_id)
        )
        if self._strand == '-':
            obj = obj.reverse_complement()
        return obj

    def extract_cds_seq(self, seqfp):
        '''
        return Sequence object
        '''
        if self._cds.is_empty():
            return None
        seq = ''.join([
            self.__extract_seq(seqfp, self._chrom, x.lower, x.upper) for x in self._cds
        ])
        obj = Sequence(
            self._id, seq, 'coding sequence, gene_id:{}'.format(self.gene_id)
        )
        if self._strand == '-':
            obj = obj.reverse_complement()
        return obj

    def extract_exon_seq(self, seqfp):
        '''
        return tuple of Sequence
        '''
        objlst = []
        exon_num = self.exon_count
        for index, x in enumerate(self._exon):
            order = index + 1 if self._strand == '+' else exon_num - index
            obj = Sequence(
                '{}_exon{}'.format(self._id, order),
                self.__extract_seq(seqfp, self._chrom, x.lower, x.upper),
                'exonic sequence, gene_id:{}'.format(self.gene_id),
            )
            if self._strand == '-':
                obj = obj.reverse_complement()
            objlst.append(obj)
        return tuple(objlst)

    def extract_intron_seq(self, seqfp):
        '''
        return tuple of Sequence
        '''
        if self.intron.is_empty():
            return tuple()
        intron_num = self.exon_count - 1    # len(self._exon)
        # assert len(self.intron) == intron_num, (self.intron, self._exon)
        objlst = []
        for index, x in enumerate(self.intron):
            # x, y = self._exon[index], self._exon[index + 1]
            order = index + 1 if self._strand == '+' else intron_num - index
            obj = Sequence(
                '{}_intron{}'.format(self._id, order),
                self.__extract_seq(seqfp, self._chrom, x.lower, x.upper),
                'intronic sequence, gene_id:{}'.format(self.gene_id),
            )
            if self._strand == '-':
                obj = obj.reverse_complement()
            objlst.append(obj)
        return tuple(objlst)

    def extract_utr5_seq(self, seqfp):
        '''
        return Sequence object
        '''
        if self._utr5.is_empty():
            return None
        seq = ''.join([
            self.__extract_seq(seqfp, self._chrom, x.lower, x.upper) for x in self._utr5
        ])
        obj = Sequence(
            '{}_utr5'.format(self._id),
            seq,
            '5\' untranslated region sequence, gene_id:{}'.format(self.gene_id)
        )
        if self._strand == '-':
            obj = obj.reverse_complement()
        return obj

    def extract_utr3_seq(self, seqfp):
        '''
        return Sequence object
        '''
        if self._utr3.is_empty():
            return None
        seq = ''.join([
            self.__extract_seq(seqfp, self._chrom, x.lower, x.upper) for x in self._utr3
        ])
        obj = Sequence(
            '{}_utr3'.format(self._id),
            seq,
            '3\' untranslated region sequence, gene_id:{}'.format(self.gene_id)
        )
        if self._strand == '-':
            obj = obj.reverse_complement()
        return obj

    def extract_utr_seq(self, seqfp):
        '''
        return tuple of Sequence
        '''
        utr5 = self.extract_utr5_seq(seqfp)
        utr3 = self.extract_utr3_seq(seqfp)
        return tuple([x for x in (utr5, utr3) if x is not None])


class GTFReader(Files):
    '''
    File parser for GTF/GFF file of gene annotation

    parameter
    ---------
    gtf:            string, file of gtf format
    flag_stream:    bool, parse file style, stream style use less memery
    '''
    __slots__ = (
        '_gene_feature', '_transcript_feature', '_utr_feature',
        '_keep_feature', '_skip_feature', '_novel_feature', '_flag_stream',
        '_drop_attr', '_strict'
    )
    def __init__(self, gtf, flag_stream=True, strict=True):
        self._gene_feature = {
            'gene', 'ncRNA_gene', 'pseudogene',
        }
        self._transcript_feature = {
            'mRNA', 'transcript',
            'tRNA', 'rRNA', 'snRNA', 'snoRNA', 'pre_miRNA',
            'pseudogenic_transcript', 'lnc_RNA', 'SRP_RNA', 'RNase_MRP_RNA'
        }
        self._utr_feature = {
            'five_prime_UTR', 'five_prime_utr', '5UTR',
            'three_prime_UTR', 'three_prime_utr', '3UTR', 'UTR',
        }
        self._keep_feature = {
            'mRNA', 'transcript', 'exon', 'CDS',
            'tRNA', 'rRNA', 'snRNA', 'snoRNA', 'pre_miRNA',
            'pseudogenic_transcript', 'lnc_RNA', 'SRP_RNA', 'RNase_MRP_RNA',
            'five_prime_UTR', 'five_prime_utr', '5UTR',
            'three_prime_UTR', 'three_prime_utr', '3UTR', 'UTR',
        }
        self._skip_feature = {
            'chromosome', 'biological_region', 'Selenocysteine',
            'start_codon', 'stop_codon',
            'gene', 'ncRNA_gene', 'pseudogene',    # skip gene feature
        }
        self._novel_feature = []
        self._flag_stream = flag_stream
        self._strict = strict
        self._drop_attr = ('ID', 'Parent', 'Name')
        Files.__init__(self, gtf)

    def __str__(self):
        return "<GTFReader object: \"{}\">".format(self._fos)

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        Iterator = self.__stream_mode() if self._flag_stream else self.__store_mode()
        for trans in Iterator:
            _transid, _chro, _interval, _strand, _exon, _cds, _utr, _attr = trans
            _attr = {
                    x: y for x, y in _attr.items() if (x not in self._drop_attr) and y
                }
            yield Transcript(
                _transid, _chro, _interval, _strand,
                exon=_exon, cds=_cds, utr=_utr, infor=_attr, strict=self._strict
            )

    def __stream_mode(self):
        '''
        '''
        t_id, t_exon, t_cds, t_utr, t_info = None, Interval(), Interval(), Interval(), {}
        for line in self.__read_gtf():
            chro, _, feature, start, end, _, strand, _, attr, line_id = line

            if t_id and ((t_id != line_id) or ((t_id == line_id) and (t_chro != chro))):
                yield (t_id, t_chro, t_interval, t_strand, t_exon, t_cds, t_utr, t_info)
                t_exon, t_cds, t_utr, t_info = Interval(), Interval(), Interval(), {}

            if feature in self._transcript_feature:
                t_id, t_chro, t_interval, t_strand = line_id, chro, intervals(start, end), strand
                try:    # gtf format
                    t_info['gene_id'] = attr.pop('gene_id')
                    t_info['gene_name'] = attr.get('gene_name', '')
                except KeyError:
                    t_info['gene_id'] = attr.pop('Parent') if ('Parent' in attr) else attr['ID']
                    t_info['transcript_name'] = attr.get('Name', '')
                t_info.update(attr)
            elif feature == 'exon':
                t_exon = t_exon | intervals(start, end)
            elif feature == 'CDS':
                t_cds = t_cds | intervals(start, end)
            elif feature in self._utr_feature:
                t_utr = t_utr | intervals(start, end)
        yield (t_id, t_chro, t_interval, t_strand, t_exon, t_cds, t_utr, t_info)

    def __store_mode(self):
        '''
        parse structure data of GTF/GFF file as python object to memory,
        store mode methods is Applicable to unsorted files,
        the same transcript annotation line is no longer in the same block
        '''
        logger.info('Start Read gff file to python object...')

        annot = OrderedDict()
        for line in self.__read_gtf():
            chro, _, feature, start, end, _, strand, _, attr, t_id = line

            if feature in self._transcript_feature:
                try:
                    geneid = attr.pop('gene_id')
                except KeyError:
                    geneid = attr.pop('Parent')
                t_info = {'gene_id': geneid}
                t_info.update(attr)
                annot.setdefault((chro, strand, t_id), {})['summary'] = [
                    intervals(start, end), t_info
                ]
            elif feature == 'exon':
                annot.setdefault((chro, strand, t_id), {}).setdefault(
                    'exon', Interval()
                ).append(intervals(start, end))
            elif feature == 'CDS':
                annot.setdefault((chro, strand, t_id), {}).setdefault(
                    'cds', Interval()
            ).append(intervals(start, end))
            elif feature in self._utr_feature:
                annot.setdefault((chro, strand, t_id), {}).setdefault(
                    'utr', Interval()
                ).append(intervals(start, end))

        logger.info('Done of parse gff, sort transcript list...')
        for uniqx in annot:
            chro, strand, transid = uniqx
            t_interval, t_info = annot[uniqx]['summary']
            t_exon = annot[uniqx].get('exon', Interval())
            t_cds = annot[uniqx].get('cds', Interval())
            t_utr = annot[uniqx].get('utr', Interval())
            yield (transid, chro, t_interval, strand, t_exon, t_cds, t_utr, t_info)

    def __read_gtf(self):
        logger.info((
            'Skip Known Annotation Feature: \n    ({})'.format(', '.join(self._skip_feature))
        ))
        for line in Files.__iter__(self):
            if line.startswith('#') or not line.strip():
                continue

            chro, source, feature, start, end, score, strand, frame, attr = line.strip().split('\t')
            if (feature not in self._keep_feature):
                if (feature not in self._skip_feature) and (feature not in self._novel_feature):
                    logger.warning('skip novel annotation feature: {}'.format(feature))
                    self._novel_feature.append(feature)
                continue
            start, end = int(start) - 1, int(end)
            if '=' in attr:
                attr = [x.strip().partition('=') for x in attr.split(';') if x.strip()]
                attr = {x[0]: x[2].strip('\'\t\n\r"') for x in attr}
                line_id = attr['ID'] if feature in self._transcript_feature else attr['Parent']
            else:
                attr = [x.strip().partition(' ') for x in attr.split(';') if x.strip()]
                attr = {x[0]: x[2].strip('\'\t\n\r"') for x in attr}
                line_id = attr['transcript_id']
            yield chro, source, feature, start, end, score, strand, frame, attr, line_id


class RefSeqGFFReader(Files):
    '''
    File parser for GFF file of gene annotation from NCBI RefSeq or Genome database

    parameter
    ---------
    gtf:     string
    chrom:   string, the table of chrom symbol convert
    strict:  bool, whether to strictly verify the gene/exon interval
    '''
    __slots__ = ('_strict', '_chrom')
    def __init__(self, gtf, chrom=None, strict=False):
        Files.__init__(self, gtf)
        self._strict = strict
        self._chrom = self.__convert_chrom_id(chrom) if chrom else {}

    def __str__(self):
        return "<RefSeqGFFReader object: \"{}\">".format(self._fos)

    def __repr__(self):
        return self.__str__()

    def __convert_chrom_id(self, tab):
        with open(tab) as f:
            tab = [i.strip().split()[:2] for i in f if not i.startswith('#')]
        return {i[0]: i[1] for i in tab}

    def __iter__(self):
        for uniqx in self.__parse_refseq_gff():
            transid_, chro, interval_, strand, exon, cds, attr = uniqx
            yield Transcript(
                transid_, self._chrom.get(chro, chro), interval_, strand,
                exon=exon, cds=cds, infor=attr, strict=self._strict
            )

    def __parse_refseq_gff(self):
        RegulateRegion = {
            'DNAseI_hypersensitive_site', 'enhancer', 'enhancer_blocking_element',
            'insulator', 'promoter', 'protein_binding_site',
            'replication_regulatory_region', 'transcriptional_cis_regulatory_region',
        }
        logger.info((
            'skip annotation feature as Regulating region: \n'
            '    ({})'.format(', '.join(RegulateRegion))
        ))

        MotifRegion = {
            'centromere', 'direct_repeat', 'microsatellite', 'tandem_repeat',
            'mobile_genetic_element', 'nucleotide_motif', 'repeat_instability_region',
        }
        logger.info((
            'skip annotation feature as Motif sequence region: \n'
            '    ({})'.format(', '.join(MotifRegion))
        ))

        UnknownRegion = {
            'cDNA_match', 'repeat_region', 'D_loop', 'match', 'region',
            'origin_of_replication', 'sequence_feature', 'biological_region',    # Igl
            'sequence_alteration', 'CAGE_cluster', 'meiotic_recombination_region',
            'mitotic_recombination_region',
        }
        logger.info((
            'skip annotation feature of unknown: \n'
            '    ({})'.format(', '.join(UnknownRegion))
        ))
        SkipRegion = UnknownRegion | MotifRegion | RegulateRegion

        GeneFeature = {'gene', 'pseudogene'}
        TransFeature = {
            'mRNA', 'lnc_RNA', 'ncRNA', 'antisense_RNA', 'transcript',  # non coding
            'RNase_MRP_RNA', 'RNase_P_RNA', 'Y_RNA', 'tRNA', 'rRNA', 'snoRNA', 'snRNA',
            'miRNA', 'primary_transcript',  # miRNA precursor seq
            'C_gene_segment', 'D_gene_segment', 'V_gene_segment',
            'J_gene_segment', 'SRP_RNA', 'telomerase_RNA',
            'vault_RNA', 'guide_RNA', 'scRNA',
        }
        logger.info((
            'the following annotation feature as transcript: \n'
            '    ({})'.format(', '.join(TransFeature))
        ))
        KeepFeature = GeneFeature | TransFeature | {'CDS', 'exon'}
        NovelFeature = []

        logger.info('Start Read gff file...')

        GeneAnnot, TransAnnot = OrderedDict(), dict()
        for line in Files.__iter__(self):
            if line.startswith('#') or not line.strip():
                continue
            chro, _, feature, start, end, _, strand, _, attr = line.strip().split('\t')
            if feature not in KeepFeature:
                if (feature not in SkipRegion) and (feature not in NovelFeature):
                    logger.warning('skip novel annotation feature: {}'.format(feature))
                    NovelFeature.append(feature)
                continue
            start, end = int(start) - 1, int(end)
            attr = [x.strip().partition('=') for x in attr.split(';') if x.strip()]
            attr = {x[0]: x[2].strip('\'\t\n\r"') for x in attr}
            # Dbxref = [i.partition(':') for i in attr.pop('Dbxref', '').split(',')]
            # Dbxref = {i[0]: i[2].strip('\'\t\n\r"') for i in Dbxref}
            # attr.update(Dbxref)

            biotype = attr.pop('gbkey', None)
            if feature in GeneFeature:
                gline_id = attr.pop('ID')
                attr['gene_name'] = attr.pop('gene')
                attr['gene_type'] = attr.pop('gene_biotype')
                GeneAnnot.setdefault((chro, strand, gline_id), {})['summary'] = [
                    intervals(start, end), attr
                ]
            elif feature in TransFeature:
                g_id = attr.pop('Parent')
                t_id = attr.pop('ID')
                if feature in ('miRNA', 'tRNA'):
                    attr['transcript_id'] = attr.pop('product')
                    attr['transcript_name'] = attr.pop('gene')
                    attr['transcript_type'] = feature
                else:
                    attr['transcript_id'] = attr.pop('transcript_id', t_id)
                    attr['transcript_name'] = attr.pop('Name', t_id)
                    attr['transcript_type'] = 'protein_coding' if biotype == 'mRNA' else biotype
                if not GeneAnnot.get((chro, strand, g_id)):
                    GeneAnnot.setdefault((chro, strand, g_id), {})['summary'] = [
                        intervals(start, end), attr
                    ]
                GeneAnnot[(chro, strand, g_id)].setdefault('translst', []).append(t_id)
                TransAnnot.setdefault(t_id, {})['summary'] = [intervals(start, end), attr]
            elif feature == 'exon':
                t_id = attr.pop('Parent')
                TransAnnot.setdefault(t_id, {}).setdefault('exon', Interval()).append(intervals(start, end))
            elif feature == 'CDS':
                t_id = attr.pop('Parent')
                TransAnnot.setdefault(t_id, {}).setdefault('cds', Interval()).append(intervals(start, end))
        logger.info('Done of Read gff, parse transcript structure ...')

        drop_attr = (
            'ID', 'Parent', 'Name', 'gene', 'gbkey', 'start_range',
            'pseudo', 'Note', 'description', 'model_evidence', 'standard_name'
        )
        index = 0
        for uniqx in GeneAnnot:
            chro, strand, geneid = uniqx
            g_interval, g_attr = GeneAnnot[uniqx]['summary']

            translst = GeneAnnot[uniqx].get('translst', [])
            flag_pseudo = False
            if not translst:
                translst.append(geneid)
                TransAnnot.setdefault(geneid, {})['summary'] = [g_interval, g_attr]
                TransAnnot[geneid]['exon'] = g_interval
                flag_pseudo = True

            for transid in translst:
                if index % 5000 == 0 and index != 0:
                    logger.info('Already processed transcript NO. : {}'.format(index))
                index += 1
                try:
                    t_interval, t_attr = TransAnnot[transid]['summary']
                except:
                    raise FileFormatError(
                        'missing transcript/gene feature in gff, tracking id: {}'.format(transid)
                    )
                t_exon = TransAnnot[transid].get('exon', Interval())
                t_cds = TransAnnot[transid].get('cds', Interval())
                assert not t_interval.is_empty()
                if not flag_pseudo:
                    t_attr['gene_id'] = geneid
                    t_attr['gene_name'] = g_attr.get('gene_name', None)
                    t_attr['gene_type'] = g_attr.get('gene_type', None)

                transid_ = t_attr.get('transcript_id', transid)
                t_attr = {
                    x: y for x, y in t_attr.items() if (x not in drop_attr) and y
                }
                yield (transid_, chro, t_interval, strand, t_exon, t_cds, t_attr)


class FileFormatError(Exception):
    def __init__(self, info):
        self.message = info


class BedReader(Files):
    '''
    '''
    __slots__ = ('_strict')
    def __init__(self, bed, strict=True):
        Files.__init__(self, bed)
        self._strict = strict

    def __str__(self):
        return "<BedReader object: \"{}\">".format(self._fos)

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        for line in Files.__iter__(self):
            line = line.strip().split('\t')
            if len(line) != 12:
                raise FileFormatError('uncomplete Bed formart, not is 12 cols')
            chro, start, end, name, _, strand, cstart, cend, _, _, elen, estart = line
            start, end, cstart, cend = map(int, [start, end, cstart, cend])
            elen = [int(x) for x in elen.split(',') if x]
            estart = [int(x) + start for x in estart.split(',') if x]
            if not len(estart)==len(elen):
                raise StructureInforError('unmatched count of exon at line: {}'.format(line))

            eend = [estart[x] + elen[x] for x, y in enumerate(elen)]
            _Exon = intervals(estart, eend)
            if cstart == cend:
                _Cds = Interval()
                attri = {'transcript_type': 'noncoding'}
            else:
                _Cds = intervals(cstart, cend) & _Exon
                attri = {'transcript_type': 'protein_coding'}
            yield Transcript(
                name, chro, intervals(start, end), strand,
                exon=_Exon, cds=_Cds, infor=attri, strict=self._strict
            )


class genePredReader(Files):
    '''
    '''
    __slots__ = ('_strict')
    def __init__(self, refFlat, strict=True):
        Files.__init__(self, refFlat)
        self._strict = strict

    def __str__(self):
        return "<genePredReader object: \"{}\">".format(self._fos)

    def __repr__(self):
        return self.__str__()

    # def __enter__(self):
    #     return self

    # def __exit__(self, rtype, value, trace):
    #     logger.debug(rtype)
    #     logger.debug(value)
    #     logger.debug(trace)

    def __iter__(self):
        for line in Files.__iter__(self):
            line = line.strip().split('\t')
            if len(line) not in (10, 11) :
                raise FileFormatError('uncomplete genePred formart, not is 10 or 11 cols')
            try:
                Tname, chro, strand, start, end, cstart, cend, ecount, estart, eend = line
                Gname = None
            except:
                Gname, Tname, chro, strand, start, end, cstart, cend, ecount, estart, eend = line
            start, end, cstart, cend, ecount = map(int, [start, end, cstart, cend, ecount])
            estart = [int(x) for x in estart.split(',') if x]
            eend = [int(x) for x in eend.split(',') if x]
            if not len(estart)==len(eend)==ecount:
                raise StructureInforError('unmatched count of exon at line: {}'.format(line))

            _Exon = intervals(estart, eend)
            if cstart == cend:
                _Cds = Interval()
                attri = {'transcript_type': 'noncoding'}
            else:
                _Cds = intervals(cstart, cend) & _Exon
                attri = {'transcript_type': 'protein_coding'}
            if Gname:
                attri['gene_id'] = Gname
                attri['gene_name'] = Gname
            yield Transcript(
                Tname, chro, intervals(start, end), strand,
                exon=_Exon, cds=_Cds, infor=attri, strict=self._strict
            )