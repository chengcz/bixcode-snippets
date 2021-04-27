#!/usr/bin/env python

import os
import sys
import requests
import pandas as pd
from pandas.io.json import json_normalize


rooturl = 'https://gnomad.broadinstitute.org/api'


def _region_vars_gnomADapi(chrom, start, end, timeout=10):

    query = """{{
  region(chrom: "{chrom}", start: {start}, stop: {end}) {{
    variants(dataset: gnomad_r2_1) {{
      variant_id: variantId
      chrom
      pos
      ref
      alt
      rsid
      reference_genome
      flags
      genome {{
        genome_ac: ac
        genome_an: an
        genome_ac_hemi: ac_hemi
        genome_ac_hom: ac_hom
        filters: filters
        # populations: populations {{
        #   id
        #   ac
        #   an
        #   ac_hom
        #   ac_hemi
        # }}
      }}
      exome {{
        exome_ac: ac
        exome_an: an
        exome_ac_hemi: ac_hemi
        exome_ac_hom: ac_hom
        filters: filters
        # populations: populations {{
        #   id
        #   ac
        #   an
        #   ac_hom
        #   ac_hemi
        # }}
      }}
    }}
  }}
}}"""
    res = requests.post(rooturl, data={'query': query.format(**locals())}, timeout=timeout)
    if res.status_code == 200:
        res = res.json()['data']['variant']
        if res is None:
            return {'variant_region': '{}:{}-{}'.format(chrom, start, end)}
        else:
            dat = json_normalize(res)
            dat = dat.fillna('.')
            return res
    print('status code is not successful, {}:{}-{}'.format(chrom, start, end))
    return None


def _variant_gnomADapi(varid, dataset, timeout=10):
    """
    varid:     chrom_pos_ref_alt
    dataset:   gnomad_r2_1
    """

    query = """{{
  variant(dataset: {dataset}, variantId: "{varid}") {{
    variant_id: variantId
    chrom
    pos
    ref
    alt
    rsid
    reference_genome
    flags
    colocatedVariants
    # sortedTranscriptConsequences {{
    #   canonical
    #   consequence_terms
    #   gene_id
    #   gene_symbol
    #   gene_version
    #   hgvs
    #   hgvsc
    #   hgvsp
    #   lof
    #   lof_flags
    #   lof_filter
    #   major_consequence
    #   polyphen_prediction
    #   sift_prediction
    #   transcript_id
    #   transcript_version
    # }}
    genome {{
      genome_ac: ac
      genome_an: an
      genome_ac_hemi: ac_hemi
      genome_ac_hom: ac_hom
      filters: filters
      populations: populations {{
        id
        ac
        an
        ac_hom
        ac_hemi
      }}
      # qualityMetrics {{
      #   siteQualityMetrics {{
      #     BaseQRankSum
      #     ClippingRankSum
      #     DP
      #     FS
      #     InbreedingCoeff
      #     MQ
      #     MQRankSum
      #     pab_max
      #     QD
      #     ReadPosRankSum
      #     RF
      #     SiteQuality
      #     SOR
      #     VQSLOD
      #   }}
      # }}
    }}
    exome {{
      exome_ac: ac
      exome_an: an
      exome_ac_hemi: ac_hemi
      exome_ac_hom: ac_hom
      filters: filters
      populations: populations {{
        id
        ac
        an
        ac_hom
        ac_hemi
      }}
      # qualityMetrics {{
      #   siteQualityMetrics {{
      #     BaseQRankSum
      #     ClippingRankSum
      #     DP
      #     FS
      #     InbreedingCoeff
      #     MQ
      #     MQRankSum
      #     pab_max
      #     QD
      #     ReadPosRankSum
      #     RF
      #     SiteQuality
      #     SOR
      #     VQSLOD
      #   }}
      # }}
    }}
  }}
}}"""
    res = requests.post(rooturl, data={'query': query.format(**locals())}, timeout=timeout)
    if res.status_code == 200:
        res = res.json()['data']['variant']
        if res is None:
            return {'variant_id': varid}
        else:
            return res
    print('status code is not successful, {}'.format(varid))
    return None

