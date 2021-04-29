
import gzip
import pandas as pd
from Bio.SeqUtils import seq1


class EffSelector(object):
    """
    parameter
    ---------
    Gene2Cano:    dict, { hgnc_id: canonical transcript, ... }
    """
    __slots__ = ('__predefined', '__distance_updown')
    def __init__(self, Gene2Cano):
        assert isinstance(Gene2Cano, dict, distance=3000)
        self.__predefined = Gene2Cano.copy()
        self.__distance_updown = distance

    def __str__(self):
        return "EffSelector object"

    def __repr__(self):
        return self.__str__()

    def annot(self, eff, fieldname, flagTab=True):
        """
        parameter
        ---------
        eff:          string, Eff infor from snpEff annot col
        fieldnames:   list,

        Eff fieldnames
        --------------
        SYMBOL, Feature, BIOTYPE, SOURCE, # REFSEQ_MATCH

        logical
        -------
        1. use "latest version" of transcript, delete older
        2. transcript type is "protein_coding"
        3. variant changes amino acid code
        4. gene symbol is "design gene of panel"
        # 5. random records from remain

        fieldname rename table
        ----------------------
          snpEff fieldname        Vep fieldname
          ----------------        -------------
          Annotation              Consequence
          Annotation_Impact       IMPACT
          Gene_Name               SYMBOL
          Gene_ID                 Gene
          Feature_Type            Feature_type
          Feature_ID              Feature
          Transcript_BioType      BIOTYPE
          Rank                    EXON/INTRON
          HGVS.c                  HGVSc
          HGVS.p                  HGVSp
          cDNA.pos / cDNA.length  cDNA_position
          CDS.pos / CDS.length    CDS_position
          AA.pos / AA.length      Protein_position
          Distance                DISTANCE
        """
        dat = pd.DataFrame(
            [x.split('|') for x in eff.split(',')], columns=fieldname
        )
        rename = {
            'Annotation':        'Consequence',
            'Annotation_Impact': 'IMPACT',
            'Gene_Name':         'SYMBOL',
            'Gene_ID':           'Gene',
            'Feature_Type':      'Feature_type',
            'Feature_ID':        'Feature',
            'Transcript_BioType': 'BIOTYPE',
            'Rank':              'EXON/INTRON',
            'HGVS.c':            'HGVSc',
            'HGVS.p':            'HGVSp',
            'cDNA.pos / cDNA.length': 'cDNA_position',
            'CDS.pos / CDS.length': 'CDS_position',
            'AA.pos / AA.length': 'Protein_position',
            'Distance':           'DISTANCE'
        }
        dat = dat.rename(rename, axis='columns')
        dat = self.__prefilrer_updown_stream(dat, self.__distance_updown)

        dat['transcript'] = dat['Feature'].apply(lambda x: x.partition('.')[0])
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
            '{},{}'.format(x, y) for x, y in zip(dat['tags'], codingtag)
        ]
        ####################
        ### canonical transcript
        self.__flag_canonical_trans(dat, self.__predefined)

        dat['tags'] = dat['tags'].apply(
            lambda x: '&'.join([m for m in x.split(',') if m.strip()])
        )
        dat['tags'] = dat['tags'].apply(lambda x: x if x else 'PASS')
        dat['HGVSp1'] = dat['HGVSp'].apply(self.__HGVSp3LetterTo1)

        if dat[dat['tags']=='PASS'].empty:
            self.__pickOne(dat)
        dat = self.__ToDict(dat, flagTab)
        return dat

    @staticmethod
    def __prefilrer_updown_stream(dt, distance):
        Over_inter = lambda x: True if \
            set([m.strip() for m in x.split('&')]) & {'intergenic_region', } else False
        if not all(map(Over_inter, dt['Consequence'])):
            dt = dt[~dt['Consequence'].apply(Over_inter)]
        else:    # all intergenic annotation
            dt = dt.iloc[[0,], ]
        if dt.shape[0] == 1:
            return dt

        updown = {'upstream_gene_variant', 'downstream_gene_variant'}
        Over_updown = lambda x: True if \
            set([m.strip() for m in x.split('&')]) & updown else False
        if all(map(Over_updown, dt['Consequence'])):
            dt['DISTANCE'] = dt['DISTANCE'].astype(int)
            if all(map(lambda x: x > distance, dt['DISTANCE'].abs())):
                dt = dt[dt['DISTANCE'].abs() == dt['DISTANCE'].abs().min()]
                dt.loc[:, 'Consequence'] = 'intergenic_region'
            else:
                dt = dt[dt['DISTANCE'].abs() == dt['DISTANCE'].abs().min()]
            dt = dt.iloc[[0,], ]
        else:    # remove upstream/downstream annotation
            dt = dt[~dt['Consequence'].apply(Over_updown)]
        return dt

    @staticmethod
    def __flag_older_version(refids):
        VERS = {}
        for refid in refids:
            ref, _, vers = refid.partition('.')
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
            chromosome_number_variation         HIGH          coding
            exon_loss_variant                   HIGH          coding
            frameshift_variant                  HIGH          coding
            rare_amino_acid_variant             HIGH          coding
            splice_acceptor_variant             HIGH          coding
            splice_donor_variant                HIGH          coding
            start_lost                          HIGH          coding
            stop_gained                         HIGH          coding
            stop_lost                           HIGH          coding
            transcript_ablation                 HIGH          coding
            3_prime_UTR_truncation & exon_loss  MODERATE      coding
            5_prime_UTR_truncation & exon_loss_variant  MODERATE  coding
            coding_sequence_variant             MODERATE      coding
            conservative_inframe_deletion       MODERATE      coding
            conservative_inframe_insertion      MODERATE      coding
            disruptive_inframe_deletion         MODERATE      coding
            disruptive_inframe_insertion        MODERATE      coding
            missense_variant                    MODERATE      coding
            regulatory_region_ablation          MODERATE
            splice_region_variant               MODERATE      coding
            TFBS_ablation                       MODERATE
            5_prime_UTR_premature_start_codon_gain_variant  LOW  coding
            initiator_codon_variant             LOW           coding
            splice_region_variant               LOW           coding
            start_retained                      LOW
            stop_retained_variant               LOW
            synonymous_variant                  LOW
            3_prime_UTR_variant                 MODIFIER
            5_prime_UTR_variant                 MODIFIER
            coding_sequence_variant             MODIFIER
            conserved_intergenic_variant        MODIFIER
            conserved_intron_variant            MODIFIER
            downstream_gene_variant             MODIFIER
            exon_variant                        MODIFIER
            feature_elongation                  MODIFIER
            feature_truncation                  MODIFIER
            gene_variant                        MODIFIER
            intergenic_region                   MODIFIER
            intragenic_variant                  MODIFIER
            intron_variant                      MODIFIER
            mature_miRNA_variant                MODIFIER
            miRNA                               MODIFIER
            NMD_transcript_variant              MODIFIER
            non_coding_transcript_exon_variant  MODIFIER
            non_coding_transcript_variant       MODIFIER
            regulatory_region_amplification     MODIFIER
            regulatory_region_variant           MODIFIER
            TF_binding_site_variant             MODIFIER
            TFBS_amplification                  MODIFIER
            transcript_amplification            MODIFIER
            transcript_variant                  MODIFIER
            upstream_gene_variant               MODIFIER
        """
        CodingRegion = {
            "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
            "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
            "inframe_insertion", "inframe_deletion", "missense_variant",
            "protein_altering_variant", "splice_region_variant",
            "incomplete_terminal_codon_variant", "start_retained_variant",
            "stop_retained_variant", "coding_sequence_variant",
            "NMD_transcript_variant", "synonymous_variant",
            # "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant"
        }
        GeneRegion = {
            "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
            "non_coding_transcript_exon_variant", "mature_miRNA_variant"
        }
        consequence = set([m.strip() for m in consequence.split('&')])
        if consequence & CodingRegion:
            return ''
        if consequence & GeneRegion:
            return 'GeneRegion'
        return 'UpDownStream'

    @staticmethod
    def __snpeff_annotation_order(annot):
        Order = {
            1: "chromosome_number_variation",
            2: "exon_loss_variant",
            3: "frameshift_variant",
            4: "stop_gained",
            5: "stop_lost",
            6: "start_lost",
            7: "splice_acceptor_variant",
            8: "splice_donor_variant",
            9: "rare_amino_acid_variant",
            10: "missense_variant",
            11: "disruptive_inframe_insertion",
            12: "conservative_inframe_insertion",
            13: "disruptive_inframe_deletion",
            14: "conservative_inframe_deletion",
            15: "5_prime_UTR_truncation+exon_loss_variant",
            16: "3_prime_UTR_truncation+exon_loss",
            17: "splice_branch_variant",
            18: "splice_region_variant",
            19: "stop_retained_variant",
            20: "initiator_codon_variant",
            21: "synonymous_variant",
            22: "initiator_codon_variant+non_canonical_start_codon",
            23: "stop_retained_variant",
            24: "coding_sequence_variant",
            25: "5_prime_UTR_variant",
            26: "3_prime_UTR_variant",
            27: "5_prime_UTR_premature_start_codon_gain_variant",
            28: "upstream_gene_variant",
            29: "downstream_gene_variant",
            30: "TF_binding_site_variant",
            31: "regulatory_region_variant",
            32: "miRNA",
            33: "custom",
            34: "sequence_feature",
            35: "conserved_intron_variant",
            36: "intron_variant",
            37: "intragenic_variant",
            38: "conserved_intergenic_variant",
            39: "intergenic_region",
            40: "coding_sequence_variant",
            41: "non_coding_exon_variant",
            42: "nc_transcript_variant",
            43: "gene_variant",
            44: "chromosome"
        }

    @staticmethod
    def __flag_canonical_trans(dat, predefined):
        # Gene, SYMBOL, transcript
        dat['_tmp'] = dat.apply(
            lambda x: predefined.get(
                x['Gene'], predefined.get(x['SYMBOL'], 'nan')
            ), axis=1
        )
        dat['tags'] = dat.apply(
            lambda x: x['tags'] if x['transcript'] == x['_tmp'] else (
                '{},{}'.format(
                    x['tags'], 'Alter' if (x['_tmp'] != 'nan') else 'unsetting'
                )
            ), axis=1
        )
        del dat['_tmp']

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

    @staticmethod
    def __pickOne(dat):
        idt = dat[dat['tags'].apply(lambda x: 'Older' not in x)]
        tmp = idt[idt['tags'].apply(lambda x: 'Noncoding' not in x)]
        if not tmp.empty: idt = tmp

        tmp = idt[idt['tags'].apply(lambda x: 'unsetting' not in x)]
        if not tmp.empty: idt = tmp

        tmp = idt[idt['tags'].apply(lambda x: 'UpDownStream' not in x)]
        if not tmp.empty: idt = tmp

        tmp = idt[idt['tags'].apply(lambda x: 'GeneRegion' not in x)]
        if not tmp.empty: idt = tmp

        tmp = idt[idt['tags'].apply(lambda x: 'Alter' not in x)]
        if not tmp.empty: idt = tmp

        if idt.shape[0] == 1:
            tag = dat.loc[idt.index[0], 'tags']
            dat.loc[idt.index, 'tags'] = 'PASS&{}'.format(tag)
        else:
            symbols = set(idt['SYMBOL'])
            for symbol in symbols:
                idx = idt[idt['SYMBOL']==symbol].index[0]
                dat.loc[idx, 'tags'] = 'random&{}'.format(dat.loc[idx, 'tags'])

    @staticmethod
    def __ToDict(dat, flagTab=True):
        dt = dat[dat['tags'].apply(lambda x: x.startswith(('PASS', 'random')))]
        if dt.empty:
            dt = dat.loc[[dat.index[0], ], ]
        dat = dat[dat['SYMBOL'].isin( set(dt['SYMBOL']) )]
        idx = list(dt.index).pop()
        dat = dat.T.to_dict()
        alts = ['Consequence', 'SYMBOL', 'Feature', 'HGVSc', 'HGVSp1', 'tags']
        Annot = dat.pop(idx)
        alter = []
        for x, y in dat.items():
            tmp = {m: y[m] for m in alts}
            tmp['HGVSc'] = tmp['HGVSc'].partition(':')[2]
            tmp['HGVSp1'] = tmp['HGVSp1'].partition(':')[2]
            alter.append(tmp)
        if flagTab:
            alter = ','.join(['|'.join([x[col] for col in alts]) for x in alter])
        Annot['alternative'] = alter
        return Annot


def read_vcf(vcf, output):
    fopen, fmode = (gzip.open, 'rt') if vcf.endswith('gz') else (open, 'r')
    with fopen(vcf, fmode) as fi:
        for line in fi:
            if line.startswith('##INFO=<ID=ANN,'):
                fieldnames = [x.strip() for x in line.split("'")[1].split('|')]
            elif line.startswith('#CHROM'):
                title = line.strip().split('\t')
                break
        dat = pd.read_csv(
            fi, sep='\t', dtype=str, comment=None, header=None, names=title
        )
        dat['POS'] = dat['POS'].astype(int)

    Selector = EffSelector(dict())
    container = []
    for _, dt in dat.iterrows():
        var = dt.to_dict()
        eff = [x for x in var['INFO'].split(';') if x.startswith('ANN=')].pop()[4:]
        var['INFO'] = ';'.join([
            x for x in var['INFO'].split(';') if not x.startswith('CSQ=')
        ])
        var.update(Selector.annot(eff, fieldnames))
        container.append(var)

    dat = pd.DataFrame(container)

