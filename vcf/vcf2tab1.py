#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import gzip
import argparse


class VariantSite(object):
    '''
    '''
    def __init__(self, line):
        chrom, pos, ID, ref, alt, _, _, infor = line.strip().split('\t')[:8]
        self._chrom = chrom
        self._pos = int(pos)
        self._id = ID
        self._ref = ref
        self._alt = alt
        self._meta = self._parse_infor(infor)

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
                # index = i.index('=')
                # infor[i[:index]] = i[index+1:]
            else:
                infor[i] = 'true'
        return infor

    def to_Tab_With_Annovar_Style(self, fp, tags):
        '''
        '''
        # 'CHROM', 'START', 'END', 'ID', 'REF', 'ALT'
        CHROM = self._chrom
        ID = self._id
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

        fp.write('\t'.join([str(i) for i in [CHROM, START, END, ID, REF, ALT]]))
        for tag in tags:
            fp.write('\t{}'.format(self._meta.get(tag, '--')))
        fp.write('\n')

    def to_Tab(self, fp, tags):
        '''
        '''
        CHROM = self._chrom
        ID = self._id
        POS = self._pos
        REF, ALT = self._ref, self._alt
        START, END = POS, POS+len(REF)-1
        fp.write('\t'.join([str(i) for i in [CHROM, POS, POS+len(REF)-1, ID, REF, ALT]]))
        for tag in tags:
            fp.write('\t{}'.format(self._meta.get(tag, '--')))
        fp.write('\n')


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
    parser.add_argument('-db', '--annovar_db',
                        dest='annovar_db', metavar='',
                        help='output annovar db file, separated by Tab')
    args = parser.parse_args()
    vcf = args.vcf
    tab = args.tab
    tags = args.tags
    tags = tags.split(',') if tags else tags
    annovar_db = args.annovar_db
    if not any([tab, annovar_db]):
        sys.stderr.write('error: \n    Missing output, --tab OR --annovar_db\n')
        sys.exit(1)
    if annovar_db:
        assert tags, '--tags must be set, output specifies INFO item'

    tab = open(tab, 'w') if tab else tab
    annovar_db = open(annovar_db, 'w') if annovar_db else annovar_db

    fopen = gzip.open if vcf.endswith('gz') else open
    with fopen(vcf) as f:
        all_tags = []
        for i in f:
            if not i.strip():
                continue

            if i.startswith('##INFO='):
                all_tags.append(i.split(',')[0].split('=')[-1])
            elif i.startswith('#CHROM'):
                all_tags = set(all_tags) - {'#CHROM', 'ALT', 'END', 'ID', 'REF', 'START'}
                tags = sorted([x for x in tags if x in all_tags] if tags else list(all_tags))
                if tab:
                    tab.write('#CHROM\tSTART\tEND\tID\tREF\tALT\t{}\n'.format('\t'.join(tags)))
                if annovar_db:
                    annovar_db.write('#CHROM\tSTART\tEND\tID\tREF\tALT\t{}\n'.format('\t'.join(tags)))
            elif not i.startswith('#'):
                mutant = VariantSite(i)
                if tab:
                    mutant.to_Tab(tab, tags)
                if annovar_db:
                    mutant.to_Tab_With_Annovar_Style(annovar_db, tags)
    if annovar_db:
        annovar_db.close()
    if tab:
        tab.close()


if __name__ == '__main__':
    main()
