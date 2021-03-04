#!/usr/bin/env python

import re
import gzip
import hashlib
import pandas as pd


def vcf2dat(vcf, meta=False, fill='.'):
    fopen, fmode = (gzip.open, 'rt') if vcf.endswith('gz') else (open, 'r')
    metainfo = []
    with fopen(vcf, fmode) as fi:
        for line in fi:
            if line.startswith('##INFO=') and meta:
                line = line.strip()[8:-1]
                line = [i.split('=') for i in line.split(',', line.count('=')-1)]
                metainfo.append({i[0]: i[1] for i in line})
            elif line.startswith('#CHROM'):
                title = line.strip('\n').split('\t')
                break
        dat = pd.read_csv(
            fi, sep='\t', dtype=str, comment=None, header=None, names=title
        )
    dat['POS'] = dat['POS'].astype(int)

    if meta:
        metainfo = pd.DataFrame(data=metainfo)
        metainfo = metainfo[['ID', 'Number', 'Type', 'Description']]
        metainfo['Description'] = metainfo['Description'].apply(lambda x: x.strip('"'))
        metas = list(set(metainfo['ID']))
        def _info_order(infor, order):
            infor = [x.split('=') for x in infor.split(';')]
            infor = {x: y for x, y in infor}
            infor = ';'.join([infor.get(x, fill) for x in order])
            return infor
        dt = dat['INFO'].apply(lambda x: _info_order(x, metas))
        dt = dt.str.split(';', expand=True)
        dt.columns = metas

        dat = pd.concat([dat, dt], axis=1, join='outer', sort=False)
        return dat, metainfo
    else:
        return dat
