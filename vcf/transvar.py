#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import re
import sys
import pysam
from io import open


def var2tab(inp, outp, refseq):
    '''
    fakevcf2tab.py
    '''
    title = [
        'gHGVS', 'input', 'type', 'Func.refGene', 'Gene.refGene',
        'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene'
    ]
    with open(refseq) as f:
        refseq = [i.strip() for i in f]

    with open(inp) as inf, open(outp, 'w') as outf:
        #print(type(outf), type(inf))
        header = '\t'.join(title)
        outf.write('chro\tpos\tref\talt\t{}\n'.format(header).decode('utf-8'))
        for line in inf:
            chro, pos, _, ref, alt, _, _, info = line.strip().split('\t')
            info = [i.split('=') for i in info.split(';')]
            info = {i[0]:i[1] for i in info if len(i)==2}
            # chr7    87138645    .    A    G    .    .    candidate_snv_variants=chr7:g.87138645A>T;candidate_codons=ATC,ATA;pHGVS=p.I1145I;gHGVS=chr7:g.87138645A>G;cHGVS=c.3435T>C;source=RefSeq;CSQN=Synonymous;reference_codon=ATT;input=ABCB1:p.I1145I;dbxref=GeneID:5243,HGNC:40,HPRD:01370,MIM:171050;type=snv;aliases=NP_000918;ANNOVAR_DATE=2018-04-16;Func.refGene=exonic;Gene.refGene=ABCB1;GeneDetail.refGene=.;ExonicFunc.refGene=synonymous_SNV;
            # AAChange.refGene=ABCB1:NM_001348946:exon26:c.T3435C:p.I1145I,ABCB1:NM_000927:exon27:c.T3435C:p.I1145I,ABCB1:NM_001348944:exon28:c.T3435C:p.I1145I,ABCB1:NM_001348945:exon30:c.T3645C:p.I1215I;ALLELE_END
            tmp = [i for i in info['AAChange.refGene'].split(',') if i.strip()]
            try:
                hgvs = [i for i in tmp if i.split(':')[1] in refseq]
            except:
                #print(tmp); break
                hgvs = []
            if not hgvs:
                hgvs = tmp
            info['AAChange.refGene'] = ','.join(hgvs)
            outf.write('{}\t{}\t{}\t{}\t{}\n'.format(chro, pos, ref, alt, '\t'.join([info.get(i, '.') for i in title])).decode('utf-8'))



def transvar2vcf(transvar, vcf):
    pattern = re.compile(r'(chr[XYM\d]+):g.([0-9]+)([ATGC]+)>([TCGA]+)')
    with open(transvar) as fi, open(vcf, 'w') as fo:
        fo.write('##fileformat=VCFv4.2\n')
        fo.write('##transvar canno --refversion hg19 --strictversion -l HGVSc.lst --gseq --refseq >variant.txt\n')
        fo.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for i in fi:
            if i.startswith('#'):
                continue
            #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
            # chrom, pos, hgvsc, ref, alt, _, _, info = i.strip().split('\t')
            inputl, transcript, gene, strand, coordinates, region, info = i.strip().split('\t')
            try:
                transcript, _ = transcript.split()
            except:
                # print(transcript)
                continue
            gDNA, cds, portein = coordinates.split('/')
            tmp = pattern.findall(gDNA)
            if tmp:
                chrom, pos, ref, alt = tmp.pop()
            else:
                print(i)
                continue
            assert inputl.split(':').pop() == cds, 'Error ... '
            attr = ('gene={};transcript={};strand={};hgvsc={};hgvsp={};'
                    '{}').format(gene, transcript, strand, cds, portein, info)
            tmp = [chrom, pos, '.', ref, alt, '.', '.', attr]
            fo.write('\t'.join(tmp) + '\n')


##############################
##############################

def parse_gHGVS(gHGVS):
    logger.debug(gHGVS)
    chro, start, end, descri_str = _split_pos_from_gHGVS(gHGVS)

    change = re.compile(r'([ATGC]*)>([ATGC]*)')
    search_change = change.search(descri_str)
    if search_change:
        ref = search_change.group(1)
        alt = search_change.group(2)
        return chro, start, end, ref, alt, 'snv'
    if descri_str.startswith('dup'):
        ref = descri_str[3:]
        return chro, start, end, ref, ref*2, 'dup'

    del_seq = re.compile(r'del([ATGC]*)')
    ins_seq = re.compile(r'ins([ATGC]*)')
    if 'ins' in descri_str and 'del' in descri_str:
        ref = del_seq.search(descri_str).group(1)
        alt = ins_seq.search(descri_str).group(1)
        return chro, start, end, ref, alt, 'delins'
    elif 'del' in descri_str:
        ref = del_seq.search(descri_str).group(1)
        return chro, start, end, ref, '', 'del'
    elif 'ins' in descri_str:
        alt = ins_seq.search(descri_str).group(1)
        return chro, start, start, '', alt, 'ins'
    else:
        return chro, start, end, '', '', 'None'


def _split_pos_from_gHGVS(gHGVS):
    chro, _, string = gHGVS.strip().partition(':g.')
    string = string.strip(' \n\r\t()_')
    for i, j in enumerate(string):
        if not j.isdigit():
            break
    start = string[:i]

    string = string[i:].strip(' \n\r\t()_')
    if string.isdigit():
        end = string
        change_descri = ''
    else:
        for i, j in enumerate(string):
            if not j.isdigit():
                break
        end = string[:i] if string[:i] else start
        change_descri = string[i:]
    return chro, int(start), int(end), change_descri


class mutation(object):
    '''
    '''
    def __init__(self, chro, pos, ref, alt, meta=None, ref_fp=None):
        self._chro = chro
        self._pos = pos if isinstance(pos, int) else int(pos)
        self._ref = ref.upper()
        self._alt = alt.upper()
        self._meta = meta
        self._fp = ref_fp
        self._valid_data()

    def _valid_data(self):
        flag_exit = []
        if set(self._ref) - set('ATGC'):
            logger.error('Error: illegal alpha in ref nucl seq: {}'.format(self._ref))
            flag_exit.append(True)
        if set(self._alt) - set('ATGC'):
            logger.error('Error: illegal alpha in alt nucl seq: {}'.format(self._alt))
            flag_exit.append(True)

        if self._fp:
            ref = self._fp.fetch(
                reference=self._chro,
                start=self._pos - 1,
                end=self._pos + len(self._ref) - 1
            )
            if ref.upper() != self._ref:
                infor = 'Error: genome seq({}) is not equal with ref seq({}) of variant'.format(ref, self._ref)
                logger.error(infor)
                flag_exit.append(True)
        if any(flag_exit):
            logger.error('Error: QC fail for variant site:')
            self.write_variant_to_vcf(sys.stderr)
            sys.exit(1)

    @property
    def chro(self):
        return self._chro

    @property
    def pos(self):
        return self._pos

    @property
    def ref(self):
        return self._ref

    @property
    def alt(self):
        return self._alt

    def write_variant_to_vcf(self, fp):
        '''
            #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
        '''
        if self._meta:
            meta = ';'.join(['{}={}'.format(i, self._meta[i]) for i in self._meta])
        else:
            meta = '.'
        fp.write(
            '{chro}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\n'.format(
                chro=self._chro,
                pos=self._pos,
                id='.',
                ref=self._ref,
                alt=self._alt,
                qual='.',
                filter='.',
                info=meta
            )
        )

    def write_variant_to_annovar(self, fp):
        '''
            #Chr     Start   End     Ref     Alt
        '''
        fp.write(
            '{chro}\t{start}\t{end}\t{ref}\t{alt}\n'.format(
                chro=self._chro,
                start=self._pos,
                end=self._pos + len(self._ref) - 1,
                ref=self._ref,
                alt=self._alt
            )
        )


def transvar2vcf(transvar, reference):
    fp_ref = pysam.FastaFile(reference)
    with open(transvar) as f:
        _ = f.readline()
        title = ['input', 'transcript', 'gene', 'strand', 'coordinates', 'region', 'info']
        for line in f:
            input_var, _, transvar_info = line.strip().partition('\t')
            transvar_info = [i.strip().split('\t') for i in transvar_info.split('|||')]
            hgvs_genome = [i[3].split('/')[0] for i in transvar_info]
            if len(set(hgvs_genome)) == 1:
                input_var_transvar = transvar_info[0]
            else:
                hgvs_genome = sorted(
                    Counter(hgvs_genome).items(), key=lambda x: x[1], reverse=True
                )
                if hgvs_genome[0][1] > hgvs_genome[1][1]:
                    gHGVS = hgvs_genome[0][0]
                    input_var_transvar = [i for i in transvar_info if i[3].split('/')[0]==gHGVS][0]
                else:
                    logger.warning('Error: Inconsistent genomic location of variant, Skip...\n    {}'.format(line))
                    # sys.exit(1)
                    continue
            transcript, gene, strand, coordinates, region, info = input_var_transvar
            if '././.' == coordinates:
                logger.debug('Debug: invalid input HGVS infor, Skip...')
                continue

            info = [i.split('=') for i in info.strip().split(';')]
            info = {i[0]: i[1] for i in info if len(i)==2}
            # CSQN=Missense;reference_codon=CAG;candidate_codons=AAG,AAA;
            # candidate_mnv_variants=chr4:g.89052321_89052323delCTGinsTTT;
            # dbxref=GeneID:9429,HGNC:74,MIM:603756;aliases=XP_005263413;source=RefSeq

            # example
            # gHGVS = info.get('candidate_mnv_variants', '').split(',')

            gHGVS, cHGVS, pHGVS = coordinates.split('/')
            chro, start, end, ref, alt, vtype = parse_gHGVS(gHGVS)
            if vtype == 'snv':
                assert len(ref)==(end-start+1), 'Error: {}'.format(line)
                pos = start
            elif vtype == 'dup':
                if not ref:
                    ref = fp_ref.fetch(reference=chro, start=start-1, end=end)
                assert len(ref)==(end-start+1), 'Error: {}'.format(line)
                pos = end
                alt = ref[-1:] + ref
                ref = ref[-1:]
            elif vtype == 'delins':
                if not ref:
                    ref = fp_ref.fetch(reference=chro, start=start-1, end=end)
                # assert len(ref)==(end-start+1), 'Error: {}'.format(line)
                if not len(ref)==(end-start+1):
                    logger.warning(
                        ('Error: length of ref seq is inconsistent '
                        'with location({}:{}-{}), Skip...\n    {}').format(
                            chro, start, end, line
                        )
                    )
                    continue
                pos = start
            elif vtype == 'del':
                if not ref:
                    ref = fp_ref.fetch(reference=chro, start=start-1, end=end)
                assert len(ref)==(end-start+1), 'Error: {}'.format(line)
                alt = fp_ref.fetch(reference=chro, start=start-2, end=start-1)
                ref = alt+ref
                pos = start - 1
            elif vtype == 'ins':
                pos = start
                ref = fp_ref.fetch(reference=chro, start=pos-1, end=pos)
                alt = ref + alt
            else:
                # chro, start, end, '', '', 'None'
                logger.warning('Warning: uncalssification variant, Skip...\n    {}'.format(line))
                continue
            tmp = {
                'type': vtype,
                'input': input_var,
                'gHGVS': gHGVS,
                'cHGVS': cHGVS,
                'pHGVS': pHGVS,
            }
            info = dict(info.items() + tmp.items())
            mutsite = mutation(
                chro,
                pos,
                ref,
                alt,
                info,
                fp_ref,    # pysam fasta file handle
            )
            yield mutsite
    fp_ref.close()



if __name__ == '__main__':
    try:
        _, transvar, reference, vcf = sys.argv
    except:
        sys.stdout.write('USAGE: \n')
        sys.stdout.write('    python script transvar.out reference.fa mutation.vcf\n\n')
        sys.stdout.write('    cmd of transvar:\n')
        sys.stdout.write('        transvar canno --refversion hg19 --refseq -l HGVS_cds --oneline >transvar.out\n')
        sys.exit(1)

    with open(vcf, 'w') as f:
        for site in transvar2vcf(transvar, reference):
            site.write_variant_to_vcf(f)
