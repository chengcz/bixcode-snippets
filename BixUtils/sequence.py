#!/usr/bin/env python

from .log import logger
from .file import Files, FileHandles, LinkFile, MakeDirs, FileName


class SequenceError(Exception):
    def __init__(self, info):
        self.message = info


class Sequence(object):
    '''
    '''
    __slots__ = ('_name', '_seq', '_descr', '_qual')
    __transtab = str.maketrans('ATGCNatgcn', 'TACGNtacgn')
    def __init__(self, name, seq, descr=None, qual=None):
        self._name = name
        self._seq = seq
        self._descr = descr if descr else ''
        if qual and len(qual) != len(seq):
            raise SequenceError('string length of sequence and qualstr is inconsistent')
        self._qual = qual

    def __str__(self):
        return "<Sequence: {}, {} ... >".format(self._name, self._seq[:50])

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return len(self._seq)

    @property
    def name(self):
        return self._name

    @property
    def seq(self):
        return self._seq

    @property
    def descr(self):
        return self._descr

    @property
    def qual(self):
        return self._qual

    @property
    def sizes(self):
        return sum([
            len(x) for x in [self._name, self._seq, self._descr, self._qual] if x
        ])

    def is_nucl(self):
        if set(self._seq) - set('ATGCNatgcn'):
            return False
        return True

    def reverse_complement(self, replace=False):
        rcseq = self._seq.translate(__transtab)[::-1]
        rcqual = self._qual[::-1] if self._qual else None
        if replace:
            self._seq = rcseq
            self._qual = rcqual
        else:
            return Sequence(self._name, rcseq, self._descr, rcqual)

    def write_to_fastq_file(self, fp):
        if self._descr:
            fp.write('@{} {}\n'.format(self._name, self._descr))
        else:
            fp.write('@{}\n'.format(self._name))
        fp.write('{}\n'.format(self._seq))
        fp.write('+\n')
        qual = self._qual if self._qual else 'I'*len(self._seq)
        fp.write('{}\n'.format(qual))

    def write_to_fasta_file(self, fp):
        if self._descr:
            fp.write('>{} {}\n'.format(self._name, self._descr))
        else:
            fp.write('>{}\n'.format(self._name))
        fp.write('{}\n'.format(self.__formater(self._seq)))

    @staticmethod
    def __formater(seq, length=80):
        fseq = ''
        for i in range(0, len(seq), length):
            fseq += seq[i : (i + length)] + '\n'
        return fseq[:-1]

    def write_to_tab_file(self, fp):
        fp.write('{}\t{}\n'.format(self._id, self._seq))


class FileFormatError(Exception):
    def __init__(self, info):
        self.message = info


class FastaReader(Files):
    '''
    File parser for fasta file of nucleic acid
    '''
    def __init__(self, fasta):
        Files.__init__(self, fasta)

    def __str__(self):
        return "<FastaReader object: \"{}\">".format(self._fos)

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        seq = None
        for line in Files.__iter__(self):
            if line.startswith('>'):
                if seq:
                    yield Sequence(seqid, seq, descr)
                seqid, _, descr = line[1:].partition(' ')
                seq = ''
            else:
                if seq is None:
                    raise FileFormatError("FASTA file does not start with '>'.")
                seq += line.strip()
        yield Sequence(seqid, seq, descr)


class FastqReader(object):
    '''
    File parser for fastq file of nucleic acid
    '''
    __slots__ = ('_fastq', '_pairedend')
    def __init__(self, fq1, fq2=None):
        self._fastq = (fq1, fq2) if fq2 else (fq1, )
        self._pairedend = (fq2 is not None)

    def __str__(self):
        return "<FastqReader object, ... >"

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        iter1 = Files(self._fastq[0]).__iter__()
        if self._pairedend:
            iter2 = Files(self._fastq[1]).__iter__()
        while True:
            try:
                name1, _, descr1 = next(iter1)[1:].partition(' ')
                if self._pairedend:
                    name2, _, descr2 = next(iter2)[1:].partition(' ')
            except StopIteration:
                break    # end of file
            try:
                seq1, _, qual1 = next(iter1), next(iter1), next(iter1)
                if self._pairedend:
                    seq2, _, qual2 = next(iter2), next(iter2), next(iter2)
            except StopIteration:
                raise FileFormatError('FASTQ file is incompelete.')
            if self._pairedend:
                if name1 != name2:
                    raise FileFormatError('Sequence order of paired fastq file is not match.')
                yield Sequence(name1, seq1, descr1, qual1), Sequence(name2, seq2, descr2, qual2)
            else:
                yield Sequence(name1, seq1, descr1, qual1)


def FastxReader(fp1=None, fp2=None, infm=None):
    '''
    parameter
    ---------
    infm:    (fasta, fastq)
    '''
    if not infm:
        infm = __guess_seq_file_format(fp1)
        if fp2 and not (infm == __guess_seq_file_format(fp2) == 'fastq'):
            raise Exception('Two file as input, but not is fastq format, now exit.')
    if infm not in ('fastq', 'fasta'):
        raise FileFormatError('Unknown file format, select from (fasta, fastq).')

    if infm == 'fasta':
        for s in FastaReader(fp1):
            yield s
    elif infm == 'fastq':
        for s in FastqReader(fp1):
            yield s


def __guess_seq_file_format(fp):
    if not fp:
        raise Exception('File Name is Null, please input right path.')

    if fp.lower().endswith(('fasta', 'fa', 'fa.gz', 'fasta.gz', 'fas')):
        return 'fasta'
    elif fp.lower().endswith(('fastq', 'fq', 'fq.gz', 'fastq.gz')):
        return 'fastq'
    else:
        raise Exception('Unknown file format, could not guess from file suffix. ')


def SplitFastxByPart(fp1, part, outdir, basename=None, fp2=None, infm=None):
    '''
    parameter
    ---------
    fp1, fp2:   file, input fasta/fastq file
    outdir      string, dir of output file
    part:       int,
    infm:       string, select from (fasta, fastq)
    '''
    if not infm:
        infm = __guess_seq_file_format(fp1)
        if fp2 and not (infm == __guess_seq_file_format(fp2) == 'fastq'):
            raise Exception('Two file as input, but not is fastq format, now exit.')
    if infm not in ('fastq', 'fasta'):
        raise FileFormatError('Unknown file format, select from (fasta, fastq).')
    if basename is None:
        basename = FileName(fp1)
    if not (part > 0):
        raise Exception('only support not negative integer')
    part = int(part)

    MakeDirs(outdir)
    if part == 1:
        logger.warning('one split file not\'s an efficient option.')
        # logger.warning('split file number is {}, do not split, soft links.'.format(part))
        # suffix = '.gz' if fp1.endswith('gz') else ''
        # if infm == 'fasta':
        #     LinkFile(fp1, '{}/{}.part01.fa{}'.format(outdir, basename, suffix))
        # elif infm == 'fastq':
        #     LinkFile(fp1, '{}/{}.part01.R1.fq{}'.format(outdir, basename, suffix))
        #     if fp2:
        #         LinkFile(fp2, '{}/{}.part01.R2.fq{}'.format(outdir, basename, suffix))
        # return

    fps = FileHandles()
    for order in range(part):
        if infm == 'fasta':
            fps.open('{}1'.format(order), '{}/{}.part{:0>3}.fa'.format(outdir, basename, order+1))
        elif infm == 'fastq':
            fps.open('{}1'.format(order), '{}/{}.part{:0>3}.R1.fq'.format(outdir, basename, order+1))
            if fp2:
                fps.open('{}2'.format(order), '{}/{}.part{:0>3}.R2.fq'.format(outdir, basename, order+1))
    for index, seq in enumerate(FastxReader(fp1, fp2, infm)):
        order = index % part
        if infm == 'fasta':
            seq.write_to_fasta_file(fps.getfp('{}1'.format(order)))
        elif infm == 'fastq' and not fp2:
            seq.write_to_fastq_file(fps.getfp('{}1'.format(order)))
        elif infm == 'fastq' and fp2:
            seq[0].write_to_fastq_file(fps.getfp('{}1'.format(order)))
            seq[1].write_to_fastq_file(fps.getfp('{}2'.format(order)))
    fps.closeall()


def SplitFastxBySizes(fp1, sizes, outdir, basename=None, fp2=None, infm=None):
    '''
    parameter
    ---------
    fp1, fp2:   file, input fasta/fastq file
    outdir      string, dir of output file
    sizes:      int,
    infm:       string, select from (fasta, fastq)
    '''
    if not infm:
        infm = __guess_seq_file_format(fp1)
        if fp2 and not (infm == __guess_seq_file_format(fp2) == 'fastq'):
            raise Exception('Two file as input, but not is fastq format, now exit.')
    if infm not in ('fastq', 'fasta'):
        raise FileFormatError('Unknown file format, select from (fasta, fastq).')
    if basename is None:
        basename = FileName(fp1)

    MakeDirs(outdir)
    order, totalbytes = 1, 0
    fps = FileHandles()
    for index, seq in enumerate(FastxReader(fp1, fp2, infm)):
        if (index == 0) or (totalbytes > sizes):
            if infm == 'fasta':
                fps.open('01', '{}/{}.part{:0>3}.fa'.format(outdir, basename, order))
            elif infm == 'fastq':
                fps.open('01', '{}/{}.part{:0>3}.R1.fq'.format(outdir, basename, order))
                if fp2:
                    fps.open('02', '{}/{}.part{:0>3}.R2.fq'.format(outdir, basename, order))
            order += 1
            totalbytes = 0

        if infm == 'fasta':
            totalbytes += seq.sizes
            seq.write_to_fasta_file(fps.getfp('01'))
        elif infm == 'fastq' and not fp2:
            totalbytes += seq.sizes
            seq.write_to_fastq_file(fps.getfp('01'))
        elif infm == 'fastq' and fp2:
            totalbytes += seq[0].sizes
            seq[0].write_to_fastq_file(fps.getfp('01'))
            seq[1].write_to_fastq_file(fps.getfp('02'))
    fps.closeall()
