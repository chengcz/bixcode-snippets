#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import gzip
import argparse


class Files(object):
    '''
    '''
    def __init__(self, File):
        self.fos = File

    def __iter__(self):
        if self.fos.lower().endswith((".gz", ".gzip")):
            lines = gzip.open(self.fos)
        else:
            lines = open(self.fos)
        for line in lines:
            try:
            #if isinstance(line, bytes):
                line = line.decode('utf-8')
            except:
                pass
            yield line
        lines.close()


class VariantSite(object):
    '''
        line: string
            one line of vcf file
        tags_info: list
            INFO ID
        tags_sample: list
            FORMAT ID
        samples: list
            sample list
    '''
    def __init__(self, line, tags_info=None, tags_format=None, samples=None):
        line = line.strip().split('\t')
        chrom, pos, ID, ref, alt, qual, filter, infor = line[:8]
        self._chrom = chrom
        self._pos = int(pos)
        self._id = ID
        self._ref = ref
        self._alt = alt
        self._qual = qual
        self._filter = filter
        self._info = self._parse_infor(infor)
        self._info_tags = tags_info
        self._format_tags = tags_format
        self._samples = samples
        self._sample_info = self._parse_sample(line[8], line[9:], samples)

    @property
    def Chrom(self):
        return self._chrom

    @property
    def Pos(self):
        return self._pos

    @property
    def ID(self):
        return self._id

    @property
    def Ref(self):
        return self._ref

    @property
    def Alt(self):
        return self._alt

    def _parse_infor(self, info):
        tmp = [i.strip() for i in info.split(';') if i.strip()]
        infor = {}
        for i in tmp:
            if '=' in i:
                i = i.split('=')
                infor[i[0]] = i[1]
            else:
                infor[i] = 'true'
        return infor

    def _parse_sample(self, Format, Info_lst, samples):
        '''
            Info_lst: lst
                vcf line, [format, sample_info1, sample_info2, ...]
            samples: list
                [sample1, sample2, ...]
        '''
        if not all([Info_lst, samples]):
            return None
        assert len(Info_lst) == len(samples), 'Error: illegal vcf format'

        sample_infor = {}
        Format = Format.split(':')
        for x, sample in enumerate(samples):
            tmp = Info_lst[x].split(':')
            sample_infor[sample] = {Format[m]: tmp[m] for m, n in enumerate(Format)}
        return sample_infor

    def to_Tab_With_Annovar_Style(self, fp, tags=None):
        '''
        '''
        if not tags:
            tags = self._info_tags

        # 'CHROM', 'START', 'END', 'ID', 'REF', 'ALT'
        POS, REF, ALT = self._pos, self._ref, self._alt
        if len(REF) == len(ALT) == 1 and REF != ALT:    # SNV
            START, END = POS, POS
        elif len(REF) == 1 and len(ALT) != 1 and ALT.startswith(REF):    # Insert
            START, END = POS, POS
            REF, ALT = '-', ALT[1:]
        elif len(REF) != 1 and len(ALT) == 1 and REF.startswith(ALT):    # Delete
            START, END = POS+1, POS+len(REF)-1
            REF, ALT = REF[1:], '-'
        else:    # complex variant
            START, END = POS, POS+len(REF)-1
        fp.write('\t'.join([str(i) for i in [self._chrom, START, END, self._id, REF, ALT]]))
        for tag in tags:
            fp.write('\t{}'.format(self._info.get(tag, '--')))
        fp.write('\n')

    def to_Tab(self, fp, tags=None, samples=None):
        '''
        '''
        if not tags:
            tags = self._info_tags
        if not samples:
            samples = self._samples

        POS, REF, ALT = self._pos, self._ref, self._alt
        START, END = POS, POS+len(REF)-1
        fp.write('\t'.join([str(i) for i in [self._chrom, POS, POS+len(REF)-1,
                                             self._id, REF, ALT, self._qual, self._filter]]))
        for tag in tags:
            fp.write('\t{}'.format(self._info.get(tag, '--')))
        if samples:
            for sample in self._samples:
                for col in self._format_tags:
                    fp.write('\t{}'.format(self._sample_info[sample].get(col, '--')))
        fp.write('\n')


class VCF_Reader(Files):
    '''
    '''
    def __init__(self, vcf, fp_tab=None, fp_adb=None, tags=None, samples=None):
        Files.__init__(self, vcf)
        infor_tags, format_tags, sample_lst = self._vcf_basic_infor()
        self._info_tags = infor_tags
        self._format_tags = format_tags
        self._sample_lst = sample_lst

        title = ['#CHROM', 'START', 'END', 'ID', 'REF', 'ALT']
        tags = self._info_tags if tags is None else tags
        if fp_tab:
            tmp = []
            samples = self._sample_lst if samples is None else samples
            for sample in samples:
                tmp += ['{}({})'.format(i, sample) for i in self._format_tags]
            fp_tab.write('{}\n'.format('\t'.join(title + ['Quality', 'Filter'] + tags + tmp)))
        if fp_adb:
            fp_adb.write('{}\n'.format('\t'.join(title + tags)))

    def __iter__(self):
        tags_info, tags_format, samples = [], [], []
        for line in Files.__iter__(self):
            if line.startswith('#'):
                continue
            else:
                yield VariantSite(line, self._info_tags, self._format_tags, self._sample_lst)

    def _vcf_basic_infor(self):
        tags_info, tags_format, samples = [], [], []
        for line in Files.__iter__(self):
            if line.startswith('##INFO='):
                tags_info.append(line.split(',')[0][11:])
            elif line.startswith('##FORMAT='):
                tags_format.append(line.split(',')[0][13:])
            elif line.startswith('#CHROM'):
                samples += line.strip().split()[9:]
        tags_info = tags_info if tags_info else None
        tags_format = tags_format if tags_format else None
        samples = samples if samples else None
        return tags_info, tags_format, samples


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf',
                        dest='vcf', metavar='',
                        required=True,
                        help='public VCF as input, [REQUIRED]')
    parser.add_argument('--tab',
                        dest='tab', metavar='',
                        help='output file separated by Tab')
    parser.add_argument('--tags',
                        dest='tags', metavar='',
                        help='select infor as columns of tab file, example: "RS,RSPOS,DP,VP"')
    parser.add_argument('-adb', '--annovar_db',
                        dest='annovar_db', metavar='',
                        help='output annovar db file, separated by Tab')
    parser.add_argument('-s', '--samples',
                        dest='samples', metavar='',
                        help='select sample for output')
    args = parser.parse_args()
    vcf = args.vcf
    tab = args.tab
    tags = args.tags
    tags = tags.split(',') if tags else tags
    annovar_db = args.annovar_db
    samples = args.samples
    samples = samples.split(',') if samples else samples

    if not any([tab, annovar_db]):
        sys.stderr.write('error: \n    Missing output, --tab OR --annovar_db\n')
        sys.exit(1)

    fp_tab = open(tab, 'w') if tab else tab
    fp_adb = open(annovar_db, 'w') if annovar_db else annovar_db
    for i in VCF_Reader(vcf, fp_tab, fp_adb, tags, samples):
        if tab:
            i.to_Tab(fp_tab, tags, samples)
        if annovar_db:
            i.to_Tab_With_Annovar_Style(fp_adb, tags)
    if tab:
        fp_tab.close()
    if annovar_db:
        fp_adb.close()


if __name__ == '__main__':
    main()
