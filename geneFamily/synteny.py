#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import pyGTF
import argparse
import numpy as np
import matplotlib.pyplot as plt
plt.rcdefaults()


def args_parse():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        epilog='''input file example:
        chrom file
        ----------
          chr1    10000

        tandem file
        -----------
          geneid1    geneid2 '''
    )
    parser.add_argument(
        '--gtf',
        dest='gtf', metavar='',
        required=True,
        help='gtf file for gene parser display, [required]'
    )
    parser.add_argument(
        '--chrom',
        dest='chrom', metavar='',
        required=True,
        help='chromosome length file, [required]'
    )
    parser.add_argument(
        '--tandem',
        dest='tandem', metavar='',
        help='optional, gene pair of tandem duplicated'
    )
    parser.add_argument(
        '-p', '--prefix',
        dest='prefix', metavar='',
        required=True,
        help='prefix of output file, [required]'
    )
    parser.add_argument(
        '-g', '--genelst',
        dest='genelst', metavar='',
        help='optional, sequence subset(1 cols) or rename(2 cols) file'
    )
    args = parser.parse_args()
    return args.gtf, args.chrom, args.tandem, args.prefix, args.genelst


def plot_scale_to_fig(length, x, y, step_y):
    R, SPACE = 0.002, 0.0025
    plt.plot(
        [x,x], [y,y-length*step_y],
        linestyle='-', color='black', lw=1.0, alpha=0.7
    )
    plt.plot(
        [x,x+2*R], [y,y],
        linestyle='-', color='black', lw=1.0, alpha=0.7
    )
    plt.plot(
        [x+2*R,x+2*R], [y,y-length*step_y],
        linestyle='-', color='black', lw=1.0, alpha=0.7
    )
    plt.plot(
        [x,x+2*R], [y-length*step_y,y-length*step_y],
        linestyle='-', color='black', lw=1.0, alpha=0.7
    )
    #plt.plot([x,x+2*R], [y,y], lw=0.5, color='black')
    for i,j in enumerate(range(5000000, int(length), 5000000)):
        plt.plot(
            [x,x+3*R], [y-j*step_y,y-j*step_y],
            color='black', lw=0.5, alpha=0.7
        )
    # unit = int(int(length)/10/1000000)    # Mb length per unit
    for i,j in enumerate(range(0, int(length), 10000000)):
        plt.plot(
            [x,x+4*R], [y-j*step_y,y-j*step_y],
            color='black', lw=0.5, alpha=0.7
        )
        plt.text(
            x+4*R+SPACE, y-j*step_y, '{:.0f} Mb'.format(j//1000000),
            color='black', fontsize=5, rotation=0, weight='semibold',
            family='Times New Roman', style='normal',
            horizontalalignment="left", verticalalignment="center"
        )


def plot_chrom_to_fig(chromlen, left, top, step_x, step_y):
    R, SPACE = 0.004, 0.05
    chr2start = {}
    for index, (chro, length) in enumerate(chromlen):
        x = left + index*step_x
        y = [top, top-length*step_y]
        plt.plot(
            [x, x], y,
            linestyle='-', color='black', lw=1.0, alpha=0.7
        )
        plt.plot(
            [x+2*R,x+2*R], y,
            linestyle='-', color='black', lw=1.0, alpha=0.7
        )
        __Telomere_halfcicle(x+R, y, R)
        plt.text(
            x+R, top+SPACE, chro,
            color='black', fontsize=7, rotation=0, weight='semibold',
            family='Times New Roman', style='normal',
            horizontalalignment="center", verticalalignment="center"
        )
        chr2start[chro] = x
    return chr2start


def __Telomere_halfcicle(x, y, r):
    top, bottom = y
    tx = np.arange(0, np.pi, np.pi/180)
    top_x = x + r*np.cos(tx)
    top_y = top + r*np.sin(tx)*2
    plt.plot(
        top_x, top_y, 'k-', color='black', lw=1.0, alpha=0.7
    )
    ty = np.arange(np.pi, 2*np.pi, np.pi/180)
    bottom_x = x + r*np.cos(ty)
    bottom_y = bottom + r*np.sin(ty)*2
    plt.plot(
        bottom_x, bottom_y, 'k-', color='black', lw=1.0, alpha=0.7
    )


def plot_gene_to_fig(chr_loc, geneloc, top, step):
    R = 0.004
    GeneID_coor, Gene_Coor = [], {}
    for (gene, chro, pos) in geneloc:
        x = chr_loc.get(chro)
        if x is None: continue
        y = top - pos*step
        plt.plot(
            [x,x+2*R], [y,y], color='black', lw=0.5, alpha=0.7
        )
        GeneID_coor.append([chro, gene, x, y])
        Gene_Coor[gene] = [x+2*R, y]
    __Label_Name(GeneID_coor, top)
    return Gene_Coor


def __Label_Name(GeneID_coor, top):
    WIDTH, SPACE, HIGHT = 0.0035, 0.002, 0.025
    chros = set([i[0] for i in GeneID_coor])
    for chro in chros:
        genelst = [i for i in GeneID_coor if i[0]==chro]
        genelst = sorted(genelst, key=lambda x: x[3], reverse=False)
        glnew = __Overlapping_Label_Adjust(
            [i[-1] for i in genelst], HIGHT, top
        )
        for i,j in enumerate(genelst):
            _, name, x1, y1 = j
            x2, y2 = x1-WIDTH, glnew[i]
            plt.plot(
                [x1,x2], [y1,y2], color='black', lw=0.5, alpha=0.7
            )
            plt.text(
                x2-SPACE, y2, name,
                color='black', fontsize=4, rotation=0,
                family='Times New Roman', style='normal',
                horizontalalignment="right", verticalalignment="center"
            )


def __Overlapping_Label_Adjust(lst, hight, top):
    pre = lst[0]
    tmp = [[pre], ]
    for i in lst[1:]:
        if i-pre > hight:
            pre = i
            tmp.append([pre, ])
        else:
            pre = pre+hight
            tmp[-1].append(pre)
    numStart = 0
    for clust in tmp:
        num = len(clust)
        numEND = numStart + num
        median_loc = lst[int((numStart+numEND)/2)]
        median = int(num/2)
        for x,y in enumerate(clust):
            clust[x] = median_loc-hight*(median-x)
        numStart = numEND
    tmp2 = []
    for i,j in enumerate(tmp):
        if i==0:
            pre = j
            tmp2.append(j)
            continue
        end, start = pre[-1], j[0]
        h = start-end
        if not h > hight:
            add_h = hight-h
            j = [x+add_h for x in j]
        tmp2.append(j)
        pre = j
    tmp = [j for i in tmp2 for j in i]
    if tmp[-1] > top:
        h = tmp[-1] - top
        tmp = [i-h for i in tmp]
    return tmp


def plot_duplicated_genepair(GeneCoor, genepair, colors='red'):
    for (gene1, gene2) in genepair:
        try:
            x1, y1 = GeneCoor[gene1]
            x2, y2 = GeneCoor[gene2]
        except:
            continue
        if x1 == x2: continue
        plt.plot(
            [x1, x2], [y1, y2], '--', lw=0.5, color=colors, #alpha=0.7
        )


def __parse_gene_location(gtf, genelst=None):
    location = []
    for t in pyGTF.GTFReader(gtf):
        if genelst and (t.id not in genelst) and (t.name not in genelst):
            continue
        tid = t.id
        if isinstance(genelst, dict):
            tid = genelst.get(t.id, genelst.get(t.name, t.id))
        location.append([tid, t.chrom, t.start])
    return location


def __parse_chrom_length(chrom, gtf):
    if chrom:
        with open(chrom) as fi:
            length = [i.strip().split() for i in fi if not i.startswith('#')]
        length = [[i[0], int(i[1])] for i in length]
    else:
        length = {}
        for t in pyGTF.GTFReader(gtf):
            if t.end > length.get(t.chrom, 0):
                length[t.chrom] = t.end
        length = [[x, y] for x, y in length.items()]
    return length


def __parse_tandem_duplicated(tandem, genelst):
    with open(tandem) as fi:
        duplicated = [i.strip().split()[:2] for i in fi if i.strip()]
    if isinstance(genelst, dict):
        duplicated = [
            (genelst[x], genelst[y]) for x, y in duplicated
            if (x in genelst) and (y in genelst)
        ]
    return duplicated


def synteny(gtf, chrom, tandem, prefix, genelst=None):
    if genelst:
        with open(genelst) as fi:
            genelst = [i.strip().split() for i in fi if i.strip()]
        genelst = {
            i[0]: i[1] for i in genelst
        } if len(genelst[0]) > 1 else {
            i[0]: i[0] for i in genelst
        }
    # else:
    #     genelst = {x[0]: x[0] for x in geneloc}

    geneloc = __parse_gene_location(gtf, genelst)
    chromlen = __parse_chrom_length(chrom, gtf)

    height, length = 0.88, 0.8
    step_x = 1.0/(len(chromlen)+1.5)
    step_y = length/max([i[1] for i in chromlen])
    start = step_x*0.9 + step_x

    plt.figure(figsize=(6.5, 3.5))
    root = plt.axes([0, 0, 1, 1])
    plot_scale_to_fig(max([i[1] for i in chromlen]), step_x*0.4, height, step_y)
    ChroStart = plot_chrom_to_fig(chromlen, start, height, step_x, step_y)
    GeneCoor = plot_gene_to_fig(ChroStart, geneloc, height, step_y)
    if tandem:
        duplicated = __parse_tandem_duplicated(tandem, genelst)
        plot_duplicated_genepair(GeneCoor, duplicated)
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()
    #plt.show()
    plt.savefig("{}.png".format(prefix), dpi=450)
    plt.savefig("{}.pdf".format(prefix))


if __name__ == '__main__':
    gtf, chrom, tandem, prefix, genelst = args_parse()
    synteny(gtf, chrom, tandem, prefix, genelst)