#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
from pyGTF import GTFReader
import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

__doc__ = '''
    ########################################
    <gtf>
    <genelst>:    order sensitive

    ########################################
'''

try:
    _, gtf, genelst, prefix = sys.argv
except:
    print('USAGE:\n    python script <gtf> <genelst> <prefix>')
    print(__doc__)
    sys.exit(1)


def fancybox(ax, dat, x, y, step, height, color):
    boxs = []
    for each in dat:
        start, end = each
        tmp = mpatches.FancyBboxPatch([x+start*step, y], step*(end-start), height,
                                      boxstyle=mpatches.BoxStyle("square", pad=0))
        boxs.append(tmp)
    collection = PatchCollection(boxs, facecolor=color,
                                 cmap=plt.cm.hsv, alpha=1, edgecolor="none")
    ax.add_collection(collection)


def splot_gene(ax, name, data, x, y, step_x):
    SPACE = 0.03
    height_intron, height_exon, height_cds = [0.0025, 0.009, 0.011]
    color_intron, color_exon, color_cds = ['black', 'blue', 'green']
    plt.text(x-SPACE, y, name, color='black', fontsize=9, rotation=0,
             family='Times New Roman', style='normal',
             horizontalalignment="right", verticalalignment="center")
    start, end, exon, cds, strand = data
    fancybox(ax, [[start-start, end-start],], x, y, step_x, height_intron, color_intron)
    # adjust gene length infor
    if exon:
        exon = [[i[0]-start, i[1]-start] for i in exon] if strand=='+' else [[end-i[1], end-i[0]] for i in exon]
        fancybox(ax, exon, x, y-(height_exon-height_intron)/2, step_x, height_exon, color_exon)
    if cds:
        cds = [[i[0]-start, i[1]-start] for i in cds] if strand=='+' else [[end-i[1], end-i[0]] for i in cds]
        fancybox(ax, cds, x, y-(height_cds-height_intron)/2, step_x, height_cds, color_cds)


def structure(gtf, genelst, prefix):
    with open(genelst) as f:
        genelst = [i.strip().split()[:2] for i in f if i.strip() and not i.startswith('#')]
    dat = {i.id:i for i in GTFReader(gtf) }#if i.id in genelst}
    length = max([dat[i].end - dat[i].start for i in dat])

    TOP, LEFT = 0.96, 0.15
    step_y = 0.94/len(dat)
    step_x = 0.82/length

    plt.figure(figsize=(7, 7))
    root = plt.axes([0, 0, 1, 1])
    for index, (gene, name) in enumerate(genelst):
        y = TOP - index * step_y
        t = dat[gene]
        StructureDat = (
            t.start, t.end,
            [(i.lower, i.upper) for i in t.exon],
            [(i.lower, i.upper) for i in t.cds],
            t.strand
        )
        splot_gene(root, name, StructureDat, LEFT, y, step_x)
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()
    #plt.show()
    plt.savefig("{}.png".format(prefix), dpi=450)
    plt.savefig("{}.pdf".format(prefix))


if __name__ == '__main__':
    structure(gtf, genelst, prefix)
