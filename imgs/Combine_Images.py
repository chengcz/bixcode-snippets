#/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import PIL.Image as Image


def combine_images(imgs, output, size = 1000, maxrow=5, maxcol=8):
    count = len(imgs)
    if count <= 3:
        row, col = 1, count
    else:
        col = int(round(pow(count, 0.5), 0))
        if col > maxrow: col = maxrow
        row = math.ceil(count / col)
        if row > maxcol: row = maxcol
    row, col = (row, col) if col > row else (col, row)
    # row, col = (row, col) if col < row else (col, row)

    if count > (maxcol * maxrow):
        print('drop some images(max 40)...')

    image = Image.new('RGB', (col * size, row * size))
    for y in range(row):
        for x in range(col):
            idx = col * y + x
            if idx >= count:
                image.paste(
                    Image.new('RGB', (col*size, row*size), (255, 255, 255)),
                    (x * size, y * size)
                )
                continue
            img = Image.open(imgs[idx])
            from_image = img.resize((size, size), Image.ANTIALIAS)
            image.paste(from_image, (x * size, y * size))
    image.save(output)


if len(sys.argv) < 4:
    print('USAGE: python %(prog)s <output.png> <img1, img2, ... >')
    sys.exit(1)

_, output, *imgs = sys.argv
combine_images(imgs, output)
