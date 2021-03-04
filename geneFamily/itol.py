#/usr/bin/env python
# coding: utf-8

import os
import sys
import glob
import pyGTF
import argparse
import pandas as pd
import xml.etree.ElementTree as ET


def __parse_genelst(genelst=None):
    if genelst:
        with open(genelst) as fi:
            geneid2name = [i.strip().split() for i in fi]
        geneid2name = {
            i[0]: i[1] for i in geneid2name
        } if len(geneid2name[0]) > 1 else {
            i[0]: i[0] for i in geneid2name
        }
        assert len(set(geneid2name.keys())) == len(set(geneid2name.values()))
    else:
        geneid2name = {}
    return geneid2name


def gtf2itol(gtf, itol, genelst, label='structure', color='#009688'):
    genelst = __parse_genelst(genelst)

    with open(itol, 'w') as fo:
        fo.write(
            'DATASET_DOMAINS\n'
            'SEPARATOR TAB\n'
            f'DATASET_LABEL\t{label}\n'
            f'COLOR\t{color}\n'
            'MARGIN\t20\n'
            'DATA\n'
        )
        for i in pyGTF.GTFReader(gtf):
            if genelst and (i.id not in genelst):
                continue
            name = genelst[i.id] if genelst else i.id

            start, end = i.start-1, i.end
            length = end - start
            tmp = []
            for color, x in zip(['#0067ff', '#ff9800', '#0067ff'], [i._utr5, i._cds, i._utr3]):
                if x.is_empty():
                    continue
                for y in x:
                    fstart, fend = y.lower, y.upper
                    if i.strand == '+':
                        tmp.append(['GP', fstart-start, fend-start, color])
                    elif i.strand == '-':
                        tmp.append(['GP', end-fend, end-fstart, color])
            tmp = sorted(tmp, key=lambda x: x[1])
            tmp = '\t'.join(
                ['{}|'.format('|'.join([str(x) for x in feature])) for feature in tmp]
            )
            fo.write('{}\t{}\t{}\n'.format(name, length, tmp))


def __parsexml_MEME(xml):
    fi = ET.parse(xml)
    root = fi.getroot()

    for item in root:
        if item.tag == 'training_set':
            inputseq = {
                i.attrib['id']: i.attrib for  i in item if i.tag == 'sequence'
            }
        elif item.tag == 'motifs':
            motifs = {}
            for motif in item:
                if motif.tag != 'motif':
                    continue
                attr = motif.attrib
                for key in motif:
                    if key.tag == 'contributing_sites':
                        for each in key:
                            attr.setdefault('sequence', []).append(each.attrib)
                motifs[attr['id']] = attr
    # reshape
    for motifid in motifs:
        mtfname = motifs[motifid]['name']
        mtflen = motifs[motifid]['width']
        for seq in motifs[motifid]['sequence']:
            seqid = seq['sequence_id']
            seqmtfpos = seq['position']
            inputseq[seqid].setdefault('motifs', []).append(
                (motifid, seqmtfpos, mtflen, mtfname)
            )
    return inputseq


def __features2visible(features):
    shapes = [
        'RE', 'HH', 'HV', 'EL', 'DI', 'TR', 'TL', 'PL', 'PR', 'PU', 'PD', 'OC', 'GP'
    ]
    colors = [
        "#377EB8", "#4DAF4A", "#FFD92F", "#E78AC3", "#A6D854",
        "#8DA0CB", "#FC8D62", "#66C2A5", "#E6AB02", "#66A61E",
        "#E7298A", "#7570B3", "#D95F02", "#1B9E77", "#F781BF",
        "#FFFF33", "#FF7F00", "#984EA3",
    ]
    while len(shapes) < len(features):
        shapes = shapes * 2
    while len(colors) < len(features):
        colors = colors * 2
    assert len(features) <= min(len(shapes), len(colors))
    feature = {i[0]: i[1:] for i in zip(features, shapes, colors)}
    return feature


def meme2itol(xml, itol, genelst, label='motif', color='#ff9800'):
    genelst = __parse_genelst(genelst)

    motifs = __parsexml_MEME(xml)
    features = set([y[0] for x in motifs for y in motifs[x].get('motifs', [])])
    feature = __features2visible(features)

    with open(itol, 'w') as fo:
        fo.write(
            'DATASET_DOMAINS\n'
            'SEPARATOR TAB\n'
            f'DATASET_LABEL\t{label}\n'
            f'COLOR\t{color}\n'
            'MARGIN\t20\n'
            'DATA\n'
        )
        for seqid in motifs:
            name = motifs[seqid]['name']
            if genelst and (name not in genelst):
                continue
            name = genelst[name] if genelst else name

            length = motifs[seqid]['length']
            tmp = []
            for motif in motifs[seqid].get('motifs', []):
                mid, start, mlen, mseq = motif
                shape, color = feature[mid]
                tmp.append('|'.join(
                    [shape, str(int(start)+1), str(int(start)+int(mlen)), color, mid]
                ))    # RE|3|24|#88f436|motif_4
            fo.write('{}\t{}\t{}\n'.format(name, length, '\t'.join(tmp)))


def __parse_smart_text(smart, skip_hidden = True):
    motifs = []
    with open(smart) as fi:
        for i in fi:
            if '=' not in i:
                continue
            values = i.strip().split('=')
            if i.startswith('DOMAIN'):
                motif = [values[1], ]
            elif i.startswith('START'):
                motif.append(values[1])
            elif i.startswith('END'):
                motif.append(values[1])
            elif i.startswith('STATUS'):
                motif.append(values[1])
                motifs.append(motif)
    motifs = [x for x in motifs if not x[-1].startswith('hidden')]
    return motifs


def smart2itol(fasta, smart_dir, itol, append=None, drop=None, genelst=None,
               label='motif', color='#673ab7'):
    genelst = __parse_genelst(genelst)

    lengths = {i.name: len(i.seq) for i in pyGTF.FastaReader(fasta)}
    if genelst:
        lengths = {genelst[x]: y for x, y in lengths.items() if genelst.get(x)}

    motifs_append = {}
    if append:
        with open(append) as fi:
            for i in fi:
                seqid_, motifid, mstart, mend = i.strip().split()
                motifs_append.setdefault(seqid_, []).append((motifid, mstart, mend, 'OK'))

    motifs = {}
    smarts = glob.glob(f'{smart_dir}/*_SMART_results.txt')
    for smart in smarts:
        seqid = os.path.basename(smart).split('_SMART')[0]
        if genelst and (seqid not in genelst):
            continue
        name = genelst[seqid] if genelst else seqid
        if not lengths.get(name, None):
            continue
        smart_mtf = __parse_smart_text(smart)
        motifs[name] = [x for x in smart_mtf if (x[0] not in drop)]
        if motifs_append.get(seqid):
            for x in motifs_append[seqid]:
                motifs[name].append(x)

    features = set([y[0] for x in motifs for y in motifs[x]])
    feature = __features2visible(features)

    with open(itol, 'w') as fo:
        fo.write(
            'DATASET_DOMAINS\n'
            'SEPARATOR TAB\n'
            f'DATASET_LABEL\t{label}\n'
            f'COLOR\t{color}\n'
            'MARGIN\t20\n'
            'DATA\n'
        )
        for name in motifs:
            length = lengths[name]
            tmp = []
            for motif in motifs[name]:
                mid, start, end, _ = motif
                shape, color = feature[mid]
                tmp.append('|'.join(
                    [shape, start, end, color, mid]
                ))    # RE|3|24|#88f436|motif_4
            fo.write('{}\t{}\t{}\n'.format(name, length, '\t'.join(tmp)))


def __exps_zscore(dat):
    rows = []
    for idx, dt in dat.iterrows():
        dt[dt.isna()] = 0
        if (dt < 1).all():
            dt = pd.Series(data=0, index=dt.index, name=dt.name)
        else:
            dt = (dt - dt.mean()) / dt.std()
        rows.append(dt)
    dat = pd.concat(rows, axis=1).T
    return dat


def exps2itol(exps, itol, genelst, label='heatmap', color='#e721f3'):
    genelst = __parse_genelst(genelst)

    dat = pd.read_csv(exps, sep='\t', index_col=0)
    rawnames = list(set(genelst.keys()) & set(dat.index))
    newnames = [genelst[x] for x in rawnames]
    dat = dat.loc[rawnames, ]
    dat.index = newnames

    dat = __exps_zscore(dat)

    color_low, color_med, color_high = '#05f705', '#ffffff', '#ff0000'
    value_low, value_med, value_high = -2, 0, 2
    color_border = '#000000'
    samples = "\t".join(dat.columns)
    with open(itol, 'w') as fo:
        fo.write(
            'DATASET_HEATMAP\n'
            'SEPARATOR TAB\n'
            f'DATASET_LABEL\t{label}\n'
            f'COLOR\t{color}\n'
            'MARGIN\t20\n'
            f'FIELD_LABELS\t{samples}\n'
            f'BORDER_WIDTH\t1\nBORDER_COLOR\t{color_border}\n'
            f'USER_MIN_VALUE\t{value_low}\nUSER_MAX_VALUE\t{value_high}\n'
            f'COLOR_MIN\t{color_low}\nCOLOR_MAX\t{color_high}\n'
            f'USER_MID_VALUE\t{value_med}\nUSE_MID_COLOR\t1\nCOLOR_MID\t{color_med}\n'
            f'COLOR_NAN\t{color_med}\n'
            'AUTO_LEGEND\t1\n'
            'LEGEND_TITLE\t\n'
            'LEGEND_SHAPES\t-2\t-1\t0\t1\t2\n'
            f'LEGEND_COLORS\t{color_low}\t#66b266\t{color_med}\t#ff6666\t{color_high}\n'
            'LEGEND_LABELS\t-2\t-1\t0\t1\t2\n'
            'DATA\n'
        )
        for idx, dt in dat.iterrows():
            geneid = dt.name
            genExps = '\t'.join([str(x) for x in dt.values])
            fo.write(f'{geneid}\t{genExps}\n')


def __align_output_formater(seq, length=None, width=60):
    if length is None:
        length = len(seq)
    fseq = ''
    for i in range(0, length, width):
        fseq += seq[i : (i + width)] + '\n'
    return fseq.strip()


def align2itol(fasta, itol, genelst, label='Alignment', color='#cddc39'):
    genelst = __parse_genelst(genelst)

    align = {s.name: s.seq for s in pyGTF.FastaReader(fasta)}
    length = [len(y) for x, y in align.items()].pop()
    with open(itol, 'w') as fo:
        fo.write(
            'DATASET_ALIGNMENT\n'
            'SEPARATOR TAB\n'
            f'DATASET_LABEL\t{label}\n'
            f'COLOR\t{color}\n'
            'CUSTOM_COLOR_SCHEME\tSCHEME_1\tA=#d2d0c9\tM=#d2d0c9\tI=#d2d0c9\t'
                'L=#d2d0c9\tV=#d2d0c9\tP=#746f69\tG=#746f69\tC=#746f69\tF=#d0ad16\t'
                'Y=#d0ad16\tW=#d0ad16\tS=#34acfb\tT=#34acfb\tN=#34acfb\tQ=#34acfb\t'
                'R=#34fb54\tK=#34fb54\tH=#34fb54\tD=#fb4034\tE=#fb4034\n'
            'CUSTOM_COLOR_SCHEME\tSCHEME_2\tA=#30a040\tR=#2015a5\tN=#10ffa0\t'
                'D=#6048c0\tC=#608a80\tQ=#601f00\tE=#5048c0\tG=#702048\tH=#a5a4a4\t'
                'I=#1a30f0\tL=#11a0f0\tK=#003505\tM=#00a0f0\tF=#10a0f0\tP=#0ff300\t'
                'S=#93f300\tT=#33ff00\tW=#0a30f0\tY=#25a4a4\tV=#90a3f0\n'
            'MARGIN\t20\n'
            'SIZE_FACTOR\t1\n'
            'COLOR_SCHEME\tclustal\n'
            # 'clustal','zappo','taylor','hphob','helix','strand','turn' and 'buried'
            'START_POSITION\t1\n'
            f'END_POSITION\t{length}\n'
            'DOTTED_DISPLAY\t0\n'
            'COLORED_DOTS\t1\n'
            'DISPLAY_CONSENSUS\t1\n'
            'DISPLAY_CONSERVATION\t1\n'
            'DATA\n'
        )
        for name, seq in align.items():
            fo.write(f'\n>{name}\n{__align_output_formater(seq)}\n')


def output_itol():
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()
    structure = subparser.add_parser('gtf2itol')
    structure.add_argument(
        '--gtf',
        dest='gtf', metavar='',
        required=True,
        help='gtf file for gene structure display, [required]'
    )
    structure.add_argument(
        '--itol',
        dest='itol', metavar='',
        required=True,
        help='output file for iTOL visualization, [required]'
    )
    structure.add_argument(
        '-g', '--genelst',
        dest='genelst', metavar='',
        help='optinal, sequence subset(1 cols) or rename(2 cols) file'
    )

    motif = subparser.add_parser('meme')
    motif.add_argument(
        '--xml',
        dest='xml', metavar='',
        required=True,
        help='XML file from MEME analysis, [required]'
    )
    motif.add_argument(
        '--itol',
        dest='itol', metavar='',
        required=True,
        help='output file for iTOL visualization, [required]'
    )
    motif.add_argument(
        '-g', '--genelst',
        dest='genelst', metavar='',
        help='optinal, sequence subset(1 cols) or rename(2 cols) file'
    )

    smart = subparser.add_parser('smart')
    smart.add_argument(
        '-f', '--fasta',
        dest='fasta', metavar='',
        required=True,
        help='fasta file for smart analysis, [required]'
    )
    smart.add_argument(
        '-indir', '--smart_results_dir',
        dest='smart_results_dir', metavar='',
        required=True,
        help='dir of smart analysis result, [required]'
    )
    smart.add_argument(
        '-a', '--append',
        dest='append', metavar='',
        help='add supplymentary motif, "seq_name\\tmotif_name\\tstart\\tend" '
    )
    smart.add_argument(
        '-d', '--drop',
        dest='drop', metavar='',
        default=['low_complexity_region', ],
        nargs='*',
        help='motif list for drop, not\'s visualization, default: low_complexity_region'
    )
    smart.add_argument(
        '--itol',
        dest='itol', metavar='',
        required=True,
        help='output file for iTOL visualization, [required]'
    )
    smart.add_argument(
        '-g', '--genelst',
        dest='genelst', metavar='',
        help='optinal, sequence subset(1 cols) or rename(2 cols) file'
    )

    heatmap = subparser.add_parser('heatmap')
    heatmap.add_argument(
        '--exps',
        dest='exps', metavar='',
        required=True,
        help='exps file with table text format, [required]'
    )
    heatmap.add_argument(
        '--itol',
        dest='itol', metavar='',
        required=True,
        help='output file for iTOL visualization, [required]'
    )
    heatmap.add_argument(
        '-g', '--genelst',
        dest='genelst', metavar='',
        help='optinal, sequence subset(1 cols) or rename(2 cols) file'
    )

    alignment = subparser.add_parser('align')
    alignment.add_argument(
        '-f', '--align-fasta',
        dest='multi_align_fasta', metavar='',
        required=True,
        help='multi alignment file with fasta format, [required]'
    )
    alignment.add_argument(
        '--itol',
        dest='itol', metavar='',
        required=True,
        help='output file for iTOL visualization, [required]'
    )
    alignment.add_argument(
        '-g', '--genelst',
        dest='genelst', metavar='',
        help='optinal, sequence subset(1 cols) or rename(2 cols) file'
    )
    submodule = ('meme', 'smart', 'gtf2itol', 'heatmap', 'align')
    if (len(sys.argv) == 1) or (sys.argv[1] not in submodule):
        parser.print_help()
        sys.exit(1)
    if sys.argv[1] == 'meme':
        args = parser.parse_args()
        meme2itol(
            args.xml, args.itol, args.genelst
        )
    elif sys.argv[1] == 'smart':
        args = parser.parse_args()
        smart2itol(
            args.fasta, args.smart_results_dir, args.itol,
            args.append, args.drop, args.genelst
        )
    elif sys.argv[1] == 'gtf2itol':
        args = parser.parse_args()
        gtf2itol(
            args.gtf, args.itol, args.genelst
        )
    elif sys.argv[1] == 'heatmap':
        args = parser.parse_args()
        exps2itol(
            args.exps, args.itol, args.genelst
        )
    elif sys.argv[1] == 'align':
        args = parser.parse_args()
        align2itol(
            args.multi_align_fasta, args.itol, args.genelst
        )


if __name__ == '__main__':
    output_itol()