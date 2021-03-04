#!/usr/bin/env python
# coding:utf-8

import os
import re
import sys
import logging
# from collections import Counter
# from os import open

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)-5s @ %(asctime)s:\n    %(message)s',
                    stream=sys.stderr)
logger = logging.getLogger(__name__)


class Variant_Canonical_Annot(object):
    '''
    parameter
    ---------
    canonical_table:    string,
        - human canonical transcript refseq ID download from
          HUGO Gene Nomenclature Committee (https://www.genenames.org),
          ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
        - non-human, input table contain two columns(gene\tcanonical transcript)
    '''

    __slots__ = ('_canonical_map')
    def __init__(self, table, hgnc=True, separator='.'):
        self._canonical_map = self.__parse_gene2transcript(table, hgnc, separator)

    def __parse_gene2transcript(self, table, flag_hgnc=False, sep='.'):
        fopen, fmode = (gzip.open, 'rt') if table.endswith('gz') else (open, 'r')
        pattern = re.compile(r'[ ,|]')

        with fopen(table, fmode) as fi:
            title = fi.readline().split()
            if flag_hgnc:
                idx_symbol = [x for x, y in enumerate(title) if y == 'symbol'].pop()
                idx_reseq = [x for x, y in enumerate(title) if y == 'refseq_accession'].pop()
            else:
                idx_symbol, idx_reseq = 0, 1
            canonical = {}
            fi.seek(0, 0)
            for line in fi:
                if line.startswith('#'):
                    continue
                line = line.split('\t')
                gene = line[idx_symbol].strip('"\t\n]r ')
                transcript = line[idx_reseq].strip('"\t\n]r ')
                if gene and transcript:
                    for x in re.split(pattern, transcript):
                        x, _ = self.__split_accession_ver(x, sep)
                        canonical.setdefault(gene, []).append(x)
        canonical = {x: tuple(canonical[x]) for x in canonical}
        return canonical

    def select_canonical(self, variant):
        '''

        selection logic
        ---------------
        - select annot from refseq database
        - biotype "protein_coding" is priority
        - pick canonical transcript Nomenclature of gene
            human canonical transcript refseq ID download from
            HUGO Gene Nomenclature Committee (https://www.genenames.org),
            ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
        - treat multi version of transcript, use lasted annot version
        - select annot according to Consequence Impact levels, ref url:
            https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html
        - select annot by distance, intergenic, upstream/downstream
        '''
        self._canonical_map
        pass

    def __vep_annot(self, annotstr, gene2acc, refseq=True, sep='.'):
        '''
        parameter
        ---------
        annotstr:   string, annotation infor from vep
        gene2acc:   dict, {symbol: canonical transcript id list of gene, ... }
        refseq:     bool, whether use refseq descr variant, default: true
        sep:        string, separator of transcript version
        '''
        title = [
            'Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type',
            'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp',
            'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons',
            'Existing_variation', 'DISTANCE', 'STRAND', 'FLAGS', 'SYMBOL_SOURCE',
            'HGNC_ID', 'REFSEQ_MATCH', 'SOURCE', 'GIVEN_REF', 'USED_REF', 'BAM_EDIT',
            'HGVS_OFFSET', 'HGVSg'
        ]
        annot = [x.split('|') for x in annotstr.split(',')]
        # step1
        if refseq:
            annot = [x for x in annot if x[24]=='RefSeq']
        # step2, biotype
        annot = self.__rank_molecular_biotype(annot)
        # step3, gene to predefined canonical accession
        annot = self.__predeflined_canonical(annot)
        # step4, multi version
        annot = self.__multi_transcrit_verison(annot)
        # step5, Consequence
        annot = self.__rank_Consequence(annot)
        # step6, closed distance in uptream / downstream / intergenic region
        annot = self.__distance_closed(annot)

        if len(annot) != 1:
            outANNstr = '\n'.join(['|'.join(x) for x in annot])
            raise CanonicalIDConfused(outANNstr)
        return annot

    @staticmethod
    def __snpeff_annot(annotstr):
        title = [
            'Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'Gene_ID',
            'Feature_Type', 'Feature_ID', 'Transcript_BioType', 'Rank',
            'HGVS.c', 'HGVS.p', 'cDNA.pos / cDNA.length', 'CDS.pos / CDS.length',
            'AA.pos / AA.length', 'Distance', 'ERRORS / WARNINGS / INFO',
        ]
        pass

    @staticmethod
    def __rank_molecular_biotype(annotlst, idx_biotype=7):
        '''
        rank by BIOTYPE
        ---------------
        eg: protein coding, biotype ref url:
        https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html
        '''
        feature = {
            '3prime_overlapping_ncrna', 'IG_C_gene', 'TR_C_gene',
            'miRNA', 'rRNA', 'snRNA', 'snoRNA', 'tRNA', 'antisense',
            'lincRNA', 'misc_RNA', 'non_stop_decay', 'nonsense_mediated_decay',
            'polymorphic_pseudogene', 'processed_pseudogene', 'pseudogene',
            'protein_coding', 'mRNA',
            'processed_transcript', 'sense_intronic', 'sense_overlapping',
            'retained_intron', 'transcribed_processed_pseudogene',
            'transcribed_unprocessed_pseudogene', 'translated_processed_pseudogene',
            'unitary_pseudogene', 'unprocessed_pseudogene',
            # ensembl
            'Mt_rRNA', 'Mt_tRNA', 'RNase_MRP_RNA', 'antisense_RNA', 'lncRNA', 'telomerase_RNA',
        }
        annot = [x for x in annotlst if x[idx_biotype]=='protein_coding']
        return annot if annot else annotlst

    @staticmethod
    def __split_accession_ver(accession, sep='.'):
            if sep in accession:
                accessionID, _, verID = accession.rpartition(sep)
                if not verID.isdigit():
                    raise AccessionError('Version is not digit.')
            else:
                accessionID, verID = accession, 0
            return accessionID, int(verID)

    # @staticmethod
    def __predeflined_canonical(self, annotlst, gene2acc, idx_symbol=3, idx_tansID=6, sep='.'):
        annot = [
            x for x in annotlst
            if self.__split_accession_ver(x[idx_tansID], sep)[0] in gene2acc.get(x[idx_symbol], ())
        ]
        return annot if annot else annotlst

    # @staticmethod
    def __multi_transcrit_verison(self, annotlst, idx_tansID=6, sep='.'):
        '''
        '''
        _MultiVer = {}
        for x in annotlst:
            transID, verID = self.__split_accession_ver(x[idx_tansID], sep)
            _MultiVer.setdefault(transID, []).append((verID, x))
        annotlst = [
            sorted(_MultiVer[x], key=lambda x: x[0], reverse=False).pop()[1]
            for x in _MultiVer
        ]
        return annotlst

    @staticmethod
    def __rank_Consequence(annotlst, idx_consequence=1):
        Consequence = {
            'intergenic_variant': 0,
            'upstream_gene_variant': 0.5, 'downstream_gene_variant': 0.5,
            'mature_miRNA_variant': 1, 'non_coding_transcript_variant': 1,
            'non_coding_transcript_exon_variant': 1,
            '3_prime_UTR_variant': 2,
            'intron_variant': 2.5, '5_prime_UTR_variant': 2.5,
            'stop_lost': 3, 'stop_retained_variant': 3, 'start_retained_variant': 3,
            'synonymous_variant': 3, 'incomplete_terminal_codon_variant': 3,
            'frameshift_variant': 4, 'inframe_deletion': 4, 'inframe_insertion': 4,
            'missense_variant': 4, 'start_lost': 4, 'stop_gained': 4,
            'coding_sequence_variant': 4, 'protein_altering_variant': 4,
            'splice_donor_variant': 4, 'splice_region_variant': 4, 'splice_acceptor_variant': 4,
            'NMD_transcript_variant': 5, 'transcript_ablation': 5,
            ### cnv and whatever
            'transcript_amplification': 3, 'TFBS_ablation': 5,
            'TFBS_amplification': 3, 'TF_binding_site_variant': 4,
            'regulatory_region_variant': 4, 'regulatory_region_ablation': 5,
            'regulatory_region_amplification': 3,
            'feature_elongation': 3, 'feature_truncation': 4,
        }
        _maxrank = max([
            Consequence[y] for x in annotlst for y in x[idx_consequence].split('&')
        ])
        annot = [
            x for x in annot
            if (max([Consequence[y] for y in x[idx_consequence].split('&')]) == _maxrank)
        ]
        return annot

    @staticmethod
    def __distance_closed(annotlst, idx_distance):
        '''
        distance, intergenic, upstream/downstream
        '''
        _distance = lambda x: int(x) if x else 0
        closed = min([_distance(x[idx_distance]) for x in annotlst])
        annot = [x for x in annotlst if _distance(x[idx_distance]) == closed]
        return annot

