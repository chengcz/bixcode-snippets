#!/usr/bin/env python

import os
import re
import sys
from pysam import FastaFile


class TransvarObj(object):
    cols = [
        'input', 'transcript', 'gene', 'strand', 'coordinates(gDNA/cDNA/protein)',
        'region', 'info'
    ]
    def __init__(self, attrs):
        if isinstance(attrs, str):
            attrs = attrs.strip().split('\t')

        if isinstance(attrs, (list, tuple)):
            ipstr, refseq, symbol, strand, hgvs, region, info = attrs
        elif isinstance(attrs, dict):
            ipstr, refseq, symbol, strand, hgvs, region, info = [attrs[x] for x in cols]
        else:
            raise Exception('unsupport type')
        self._rawip = (ipstr, refseq, symbol, strand, hgvs, region, info)
        self.id = ipstr
        self.hgvsg, self._hgvsc, self._hgvsp = hgvs.split('/')

        fieldname = ['CSQN', 'dbxref', 'aliases', 'candidate_mnv_variants']
        self.consequence, dbxref, self.aliases, self.alternative = self._parse_info(info, fieldname)

        fieldname = ['GeneID', 'HGNC', 'MIM']
        self.geneid, self.hgncid, self.omim = self._parse_xref(dbxref, fieldname)
        self.symbol, self.refseq = symbol, refseq.split()[0]
        self.type = self._type(self.hgvsg)
        self.__flagloc = False

    def __str__(self):
        msg  = '<TransvarObj object: {}; {}>'.format(self.id, self.hgvsg)
        return msg

    def __repr__(self):
        return self.__str__()

    def is_valid(self):
        if self.hgvsg == '.':
            return False
        return True

    @staticmethod
    def _parse_info(info, field):
        info = [x.partition('=') for x in info.split(';')]
        info = {x: y if y else True for (x, _, y) in info}
        return [info.get(x, '.') for x in field]

    @staticmethod
    def _parse_xref(xref, field):
        tmp = [m.partition(':') for m in xref.split(',')]
        xref = {x: y for (x, _, y) in tmp}
        return [xref.get(x, '.') for x in field]

    @staticmethod
    def _type(hgvsg):
        if '>' in hgvsg:
            return 'snv'
        elif ('del' in hgvsg) and ('ins' in hgvsg):
            return 'delins'
        elif ('del' in hgvsg):
            return 'del'
        elif ('ins' in hgvsg):
            return 'ins'
        elif ('dup' in hgvsg):
            return 'dup'
        else:
            return None

    @property
    def hgvsc(self):
        if not self.is_valid():
            return '.'
        return '{}:{}'.format(self.refseq, self._hgvsc)

    @property
    def hgvsp(self):
        if not self.is_valid():
            return '.'
        return '{}:{}'.format(self.aliases, self._hgvsp)

    @property
    def location(self):
        if self.__flagloc:
            return self.chrom, self.pos, self.ref, self.alt

        if not self.is_valid():
            self.chrom, self.pos, self.ref, self.alt = '.', '.', '.', '.'
            return self.chrom, self.pos, self.ref, self.alt

        self.chrom, self.pos, self.ref, self.alt = self._hgvsg_to_chroms(self.hgvsg)
        self.__flagloc = True
        return self.chrom, self.pos, self.ref, self.alt

    @property
    def allLocation(self):
        if not self.__flagloc:
            self.location

        item = [(self.chrom, self.pos, self.ref, self.alt), ]
        if self.alternative != '.':
            alters = [self._hgvsg_to_chroms(x.strip()) for x in self.alternative.split(',')]
            item += alters
        return tuple(item)

    @staticmethod
    def _hgvsg_to_chroms(hgvsg):
        def nonDigitIdx(string):
            for idx, letter in enumerate(string):
                if not letter.isdigit():
                    break
            return idx

        chrom, hgvsg_suffix = hgvsg.split(':g.')
        if '>' in hgvsg_suffix:
            pos_ref, alt = hgvsg_suffix.split('>')
            idx = nonDigitIdx(pos_ref)
            pos, ref = pos_ref[:idx], pos_ref[idx:]
        elif ('del' in hgvsg_suffix) and ('ins' in hgvsg_suffix):
            idel, iins = hgvsg_suffix.index('del'), hgvsg_suffix.index('ins')
            pos, *end = hgvsg_suffix[:min(idel, iins)].split('_')
            ref = hgvsg_suffix[idel+3:iins]
            alt = hgvsg_suffix[iins+3:]
            if ref == '':    # chr1:g.10689923_10689924delinsGT
                assert len(end) in (0, 1)
                ref = 'N' * (int(end[0]) - int(pos) + 1) if end else 'N'
        elif ('del' in hgvsg_suffix):
            dix = hgvsg_suffix.index('del')
            ref = hgvsg_suffix[dix+3:]
            pos, *end = hgvsg_suffix[:dix].split('_')
            flag = False
            for idx, letter in enumerate(ref.upper()):
                if not (letter in set('ATGCN')):
                    flag = True
                    break
            if flag: ref = ref[:idx]
            if ref == '':    # chr7:g.55242467_55242481del15
                assert len(end) in (0, 1)
                ref = 'N' * (int(end[0]) - int(pos) + 1) if end else 'N'
            alt = '.'
        elif ('ins' in hgvsg_suffix):
            idx = nonDigitIdx(hgvsg_suffix)
            pos = hgvsg_suffix[:idx]
            alt = hgvsg_suffix[hgvsg_suffix.index('ins')+3:]
            ref = '.'
        elif ('dup' in hgvsg_suffix):
            dix = hgvsg_suffix.index('dup')
            alt = hgvsg_suffix[dix+3:]
            pos, *end = hgvsg_suffix[:dix].split('_')
            if alt:
                ref = alt[:1]
                alt = alt + ref
            else:
                assert len(end) in (0, 1)
                ref = 'N' * (int(end[0]) - int(pos) + 1) if end else 'N'
                alt = ref + ref
        else:
            chrom, pos, ref, alt = '.', 0, '.', '.'
        return (chrom, int(pos), ref, alt)

    def validate_Ref(self, fai):
        """
        not implemented
        """
        raise NotImplementedError('not implemented')
        if not self.__flagloc:
            self.location
        chrom, pos, ref, alt = self.chrom, self.pos, self.ref, self.alt
        if chrom.startswith('chr'):
            chrom = chrom[3:]

        # fai = pysam.FastaFile(reference)
        if ref.startswith('N'):    # ins, delins, del
            ref = fai.fetch(reference=chrom, start=pos-1, end=pos+len(ref)-1)
        if alt == '.':   # del
            pos = pos - 1
            alt = fai.fetch(reference=chrom, start=pos-1, end=pos)
            ref = alt + ref
        ref = ref.upper()
        self.pos, self.ref, self.alt = pos, ref, alt.upper()

        REF = fai.fetch(reference=chrom, start=pos-1, end=pos+len(ref)-1).upper()
        return ref == REF

    def vcf(self):
        if not self.__flagloc:
            self.location
        attrs = {
            'hgvsg': self.hgvsg,
            'hgvsc': self.hgvsc,
            'hgvsp': self.hgvsp,
            'consequence': self.consequence,
            'symbol': self.symbol,
            'refseq': self.refseq,
            'aliases': self.aliases,
            'geneid': self.geneid,
            'hgncid': self.hgncid,
            'omim': self.omim,
            'alternative': self.alternative,
        }
        return [
            self.chrom, self.pos, self.id, self.ref, self.alt, '.', '.',
            ';'.join(['{}={}'.format(x, y) for x, y in attrs.items()])
        ]


class VcfWriter(object):
    def __init__(self, fp):
        self._fp = open(fp, 'r')
        self.__header()

    def __header(self):
        header  = '##fileformat=VCFv4.1\n'
        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        self._fp.write(header)

    def write(self, content, sep='\t'):
        if isinstance(content, (list, tuple)):
            content = sep.join(map(str, content))
        elif isinstance(content, str):
            pass
        else:
            raise Exception('Unsupport Type')
        self._fp.write(content)

    def close(self):
        self._fp.close()


##############################
##### old methods
##############################
def parseHGVSg(HGVSg):
    chro, start, end, var_suffix = __extract_location(HGVSg)

    patternSnv = re.compile(r'([ATGC]*)>([ATGC]*)')
    snvmatch = patternSnv.search(var_suffix)
    if snvmatch:
        ref = snvmatch.group(1)
        alt = snvmatch.group(2)
        return chro, start, end, ref, alt, 'snv'
    if var_suffix.startswith('dup'):
        ref = var_suffix[3:]
        return chro, start, end, ref, ref*2, 'dup'

    patternDel = re.compile(r'del([ATGC]*)')
    patternIns = re.compile(r'ins([ATGC]*)')
    if 'ins' in var_suffix and 'del' in var_suffix:
        ref = patternDel.search(var_suffix).group(1)
        alt = patternIns.search(var_suffix).group(1)
        return chro, start, end, ref, alt, 'delins'
    elif 'del' in var_suffix:
        ref = patternDel.search(var_suffix).group(1)
        return chro, start, end, ref, '', 'del'
    elif 'ins' in var_suffix:
        alt = patternIns.search(var_suffix).group(1)
        return chro, start, start, '', alt, 'ins'
    else:
        return chro, start, end, '', '', 'nan'


def __extract_location(HGVSg):
    chro, _, string = HGVSg.strip().partition(':g.')
    string = string.strip(' \n\r\t()_')
    for i, j in enumerate(string):
        if not j.isdigit():
            break
    start = string[:i]

    string = string[i:].strip(' \n\r\t()_')
    if string.isdigit():
        end = string
        suffix = ''
    else:
        for i, j in enumerate(string):
            if not j.isdigit():
                break
        end = string[:i] if string[:i] else start
        suffix = string[i:]
    return chro, int(start), int(end), suffix
