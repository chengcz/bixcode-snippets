#!/usr/bin/env python

import os
import sys
import gzip
import json
from io import open
import pandas as pd
from Bio.SeqUtils import seq1


class CsqSelector(object):
    """
    parameter
    ---------
    transcript:    path to file, two columns, "symbol \\t refseq"
    hgnc:          path to hgnc.txt, download genenames.org
    refFlat:       path to refFlat.txt, download ucsc.com
    """
    __slots__ = ('__input', '__predefined', '__hgnc', '__reflen')
    def __init__(self, transcript=None, hgnc=None, refFlat=None):
        self.__predefined = self.__parse_defined(transcript) if transcript else None
        self.__hgnc = self.__parse_hgnc(hgnc) if hgnc else None
        self.__reflen = self.__parse_refFlat(refFlat) if refFlat else None
        self.__input = {'defined': transcript, 'hgnc': hgnc, 'refFlat': refFlat}

    def __str__(self):
        msg  = '='*60 + '\n'
        msg += 'CsqSelector object\n'
        msg += '  define symbol/hgncid to canonical transcript, reference:\n'
        msg += '  >> predefined transcript: {}\n'.format(self.__input['defined'])
        msg += '  >> complete hgnc:         {}\n'.format(self.__input['hgnc'])
        msg += '  >> refseq gene predict:   {}\n'.format(self.__input['refFlat'])
        msg += '='*60 + '\n'
        return msg

    def __repr__(self):
        return self.__str__()

    @staticmethod
    def __parse_defined(text):
        fopen, fmode = (gzip.open, 'rt') if text.endswith('gz') else (open, 'r')
        sym2refid = {}
        with fopen(text, fmode) as fi:
            # header = fi.readline().strip('\r\n').split('\t')
            for line in fi:
                symbol, refseq = line.strip('\r\n').split()[:2]
                # line = line.strip('\r\n').split('\t')
                # line = {x: y for x, y in zip(header, line)}
                # sym2refid[line['symbol']] = line['refseq']
                sym2refid[symbol] = refseq
        return sym2refid

    @staticmethod
    def __parse_hgnc(text):
        """
        Canonical transcript refseq ID of gene Nomenclature download from
            HUGO Gene Nomenclature Committee (https://www.genenames.org),
            ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
        """
        fopen, fmode = (gzip.open, 'rt') if text.endswith('gz') else (open, 'r')
        hgnc2refseq = {}
        with fopen(text, fmode) as fi:
            header = fi.readline().strip('\r\n').split('\t')
            for line in fi:
                line = line.strip('\r\n').split('\t')
                line = {x: (y if y else '.') for x, y in zip(header, line)}
                if line['refseq_accession'] == '.':
                    continue
                hgnc_id = line['hgnc_id'][5:]
                symbol  = line['symbol']
                refseq  = line['refseq_accession'].split('|')
                hgnc2refseq[hgnc_id] = refseq
                hgnc2refseq[symbol]  = refseq
        return hgnc2refseq

    @staticmethod
    def __parse_refFlat(text):
        fopen, fmode = (gzip.open, 'rt') if text.endswith('gz') else (open, 'r')
        # NACA2  NM_199290  chr17  -  59667783  59668580  59667893  59668541  1  59667783,  59668580,
        refseq2len = {}
        str2int = lambda x: map(int, [m for m in x.split(',') if m])
        with fopen(text, fmode) as fi:
            for line in fi:
                # symbol, rsid, chrom, strand, start, end, \
                # txstart, txend, count, exstart, exend
                line = line.strip().split('\t')
                exstart, exend = map(str2int, [line[-2], line[-1]])
                exlen = sum([end-start for start, end in zip(exstart, exend)])
                refid = line[1]
                if exlen > refseq2len.get(refid, 0):
                    refseq2len[refid] = exlen
        return refseq2len

    def annot(self, csq, fieldname):
        """
        parameter
        ---------
        csq:          string, CSQ infor from VEP annot
        fieldnames:   list,

        CSQ fieldnames
        --------------
        SYMBOL, Feature, BIOTYPE, SOURCE, # REFSEQ_MATCH

        logical
        -------
        1. whether or not transcript ID from "RefSeq"
        2. use "latest version" of transcript, delete older
        3. transcript type is "protein_coding"
        4. variant changes amino acid code, ref as
            https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html
        5. gene symbol is "design gene of panel"
        # 6. "rna match" records if refseq
        # 7. random records from remain
        """
        dat = pd.DataFrame(
            [x.split('|') for x in csq.split(',')], columns=fieldname
        )
        tmp = dat[dat['SOURCE']=='RefSeq']
        if not tmp.empty: dat = tmp
        dat['transcript'] = dat['Feature'].apply(lambda x: x.rpartition('.')[0])

        dat['tags'] = ''
        ### tag older version
        OLDER = self.__flag_older_version(dat['Feature'])
        dat['tags'] = dat.apply(
            lambda x: '{},Older'.format(x['tags']) \
                if (x['Feature'] in OLDER) else x['tags'],
            axis=1
        )
        ### protein coding
        dat['tags'] = dat.apply(
            lambda x: '{},Noncoding'.format(x['tags']) \
                if (x['BIOTYPE'] != 'protein_coding') else x['tags'],
            axis=1
        )
        ### coding region
        codingtag = dat['Consequence'].apply(self.__flag_coding_region)
        dat['tags'] = [
            ','.join([m for m in [x, y] if m.strip()])
            for x, y in zip(dat['tags'], codingtag)
        ]
        ### canonical transcript
        self.__flag_canonical_trans(dat, self.__predefined, self.__hgnc, self.__reflen)

        dat['tags'] = dat['tags'].apply(lambda x: x.strip(',') if x else 'PASS')
        dat['HGVSp1'] = dat['HGVSp'].apply(self.__HGVSp3LetterTo1)

        if dat[dat['tags']=='PASS'].empty:
            self.__pickOne(dat)
        dat = self.__ToDict(dat)
        return dat

    @staticmethod
    def __ToDict(dat):
        dt = dat[dat['tags'].apply(lambda x: x.startswith(('PASS', 'random')))]
        if dt.empty:
            dt = dat.loc[[dat.index[0], ], ]
        dat = dat[dat['SYMBOL'].isin( set(dt['SYMBOL']) )]
        idx = list(dt.index).pop()
        cols = [
            'Consequence', 'SYMBOL', 'Gene', 'HGNC_ID', 'STRAND',
            'Feature', 'BIOTYPE', 'EXON', 'INTRON',
            'HGVSg', 'HGVSc', 'HGVSp', 'HGVSp1', 'tags',
            'CDS_position', 'Codons', 'Protein_position', 'Amino_acids',
        ]
        dat = dat[cols].T.to_dict()
        alts = [
            'Consequence', 'SYMBOL', 'Gene', 'HGNC_ID', 'HGVSc', 'HGVSp1', 'tags'
        ]
        Annot = dat.pop(idx)
        Annot['alternative'] = tuple([
            {m: y[m] for m in alts} for x, y in dat.items()
        ])
        return Annot

    @staticmethod
    def __pickOne(dat):
        idt = dat[dat['tags'].apply(lambda x: 'Older' not in x)]
        tmp = idt[idt['tags'].apply(lambda x: 'Noncoding' not in x)]
        if not tmp.empty: idt = tmp

        tmp = idt[idt['tags'].apply(lambda x: 'UpDownStream' not in x)]
        if tmp.empty:
            dist = str(idt['DISTANCE'].astype(int).min())
            glst = set(idt.loc[idt['DISTANCE']==dist, 'SYMBOL'])
            idt = idt[idt['SYMBOL'].isin(glst)]
        else:
            idt = tmp

        tmp = idt[idt['tags'].apply(lambda x: 'GeneRegion' not in x)]
        if not tmp.empty: idt = tmp

        tmp = idt[idt['tags'].apply(lambda x: 'Alter' not in x)]
        if not tmp.empty: idt = tmp

        if idt.shape[0] == 1:
            tag = dat.loc[idt.index[0], 'tags']
            dat.loc[idt.index, 'tags'] = 'PASS,{}'.format(tag)
        else:
            symbols = set(idt['SYMBOL'])
            for symbol in symbols:
                idx = idt[idt['SYMBOL']==symbol].index[0]
                dat.loc[idx, 'tags'] = 'random,{}'.format(dat.loc[idx, 'tags'])
        return

    @staticmethod
    def __flag_older_version(refids):
        VERS = {}
        for refid in refids:
            ref, _, vers = refid.rpartition('.')
            VERS.setdefault(ref, []).append((
                refid, int(vers) if vers.isdigit() else 1
            ))
        delrefids = []
        for _, refids in VERS.items():
            *tmp, _ = sorted(refids, key=lambda m: m[1])
            delrefids.extend(tmp)
        return tuple([x[0] for x in delrefids])

    @staticmethod
    def __flag_coding_region(consequence):
        """
        reference
        ---------
            https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html

        Sequence Ontology term
        ----------------------
            transcript_ablation                       HIGH              coding
            splice_acceptor_variant                   HIGH              coding
            splice_donor_variant                      HIGH              coding
            stop_gained                               HIGH              coding
            frameshift_variant                        HIGH              coding
            stop_lost                                 HIGH              coding
            start_lost                                HIGH              coding
            transcript_amplification                  HIGH
            inframe_insertion                         MODERATE          coding
            inframe_deletion                          MODERATE          coding
            missense_variant                          MODERATE          coding
            protein_altering_variant                  MODERATE          coding
            splice_region_variant                     LOW               coding
            incomplete_terminal_codon_variant         LOW               coding
            start_retained_variant                    LOW               coding
            stop_retained_variant                     LOW               coding
            synonymous_variant                        LOW               @ coding
            coding_sequence_variant                   MODIFIER          coding
            mature_miRNA_variant                      MODIFIER
            5_prime_UTR_variant                       MODIFIER
            3_prime_UTR_variant                       MODIFIER
            non_coding_transcript_exon_variant        MODIFIER
            intron_variant                            MODIFIER
            NMD_transcript_variant                    MODIFIER          coding
            non_coding_transcript_variant             MODIFIER
            upstream_gene_variant                     MODIFIER
            downstream_gene_variant                   MODIFIER
            TFBS_ablation                             MODIFIER
            TFBS_amplification                        MODIFIER
            TF_binding_site_variant                   MODIFIER
            regulatory_region_ablation                MODERATE
            regulatory_region_amplification           MODIFIER
            feature_elongation                        MODIFIER
            regulatory_region_variant                 MODIFIER
            feature_truncation                        MODIFIER
            intergenic_variant                        MODIFIER
        """
        Coding = {
            "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
            "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
            "inframe_insertion", "inframe_deletion", "missense_variant",
            "protein_altering_variant", "splice_region_variant",
            "incomplete_terminal_codon_variant", "start_retained_variant",
            "stop_retained_variant", "coding_sequence_variant",
            "NMD_transcript_variant", "synonymous_variant",
            # "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant"
        }
        if set(consequence.split('&')) & Coding:
            return ''
        Gene = {
            "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
            "non_coding_transcript_exon_variant", "mature_miRNA_variant"
        }
        if set(consequence.split('&')) & Gene:
            return 'GeneRegion'
        return 'UpDownStream'

    @staticmethod
    def __flag_canonical_trans(dat, s_predefined=None, s_hgnc=None, s_reflen=None):
        flag = False
        dt = dat[['SYMBOL', 'HGNC_ID']].drop_duplicates()
        if s_predefined:
            dt['transcript'] = dt['SYMBOL'].apply( lambda x: s_predefined.get(x) )
            tmp = dt[~dt['transcript'].isna()]
            if not tmp.empty:
                transids = set(tmp['transcript'])
                if set(dat['transcript']) & transids:
                    dat['tags'] = dat.apply(
                        lambda x: '{},AlterDefine'.format(x['tags']) \
                            if (x['transcript'] not in transids) else x['tags'],
                        axis=1
                    )
                    return
        if s_hgnc:
            dt['transcript'] = dt.apply(
                lambda x: s_hgnc.get(x['HGNC_ID'], s_hgnc.get(x['SYMBOL'])), axis=1
            )
            tmp = dt[~dt['transcript'].isna()]
            if not tmp.empty:
                transids = set([y for x in tmp['transcript'] for y in x])
                final = set(dat['transcript']) & transids
                if final:
                    if len(final) > 1:
                        transids_ = set([x[0] for x in tmp['transcript']])
                        if set(dat['transcript']) & transids_:
                            transids = transids_
                    dat['tags'] = dat.apply(
                        lambda x: '{},AlterHGNC'.format(x['tags']) \
                            if (x['transcript'] not in transids) else x['tags'],
                        axis=1
                    )
                    return
        if s_reflen:
            maxlen = dat['transcript'].apply(lambda x: s_reflen.get(x, 0)).max()
            dat['tags'] = dat.apply(
                lambda x: '{},AlterLens'.format(x['tags']) \
                    if s_reflen.get(x['transcript']) != maxlen else x['tags'],
                axis=1
            )
            return
        return

    @staticmethod
    def __HGVSp3LetterTo1(hgvsp):
        Cache, flag = [], None
        lower = lambda x: 'a' <= x <= 'z'
        if (':' in hgvsp):
            pid, _, hgvsp =  hgvsp.partition(':')
            Cache.extend([pid, ':'])

        hgvsp = hgvsp.replace('%3D', '=', 1) if ('%3D' in hgvsp) else hgvsp
        for idx, x in enumerate(hgvsp):
            if flag and (idx < flag): continue
            if 'A' <= x <= 'Z':
                if all(map(lower, hgvsp[idx+1:idx+3])):
                    Cache.append(seq1(hgvsp[idx:idx+3]))
                    flag = idx + 3
                else:
                    Cache.append(x)
            else:
                Cache.append(x)
        return ''.join(Cache)


try:
    _, hgnc, refFlat, vcf, prefix = sys.argv
except:
    print('usage: python %(prog)s hgnc refflat vcf prefix')
    sys.exit(1)


selector = CsqSelector(hgnc=hgnc, refFlat=refFlat)

with open(vcf) as fi:
    lines = []
    for idx, line in enumerate(fi):

        if line.startswith('##INFO=<ID=CSQ'):
           fieldname = line.split(': ').pop().strip('">\n').split('|')
        elif line.startswith('#CHROM'):
           title = line.strip().split()
        if line.startswith('#'):
           continue

        line = {x: y for x, y in zip(title, line.strip().split('\t'))}
        infor = [x.partition('=') for x in line['INFO'].split(';')]
        infor = {x: y for (x, _, y) in infor}

        try:
            line['CSQ'] = selector.annot(infor['CSQ'], fieldname)
        except:
            line['CSQ'] = 'error'
