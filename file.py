#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import bz2
import sys
import gzip
import logging


logging.basicConfig(
    stream=sys.stderr
    level=logging.INFO,
    format='%(levelname)-5s @ %(filename)s:%(lineno)s, %(asctime)s, %(message)s',
)
logger = logging.getLogger(__name__)
# stdlog = logging.StreamHandler(sys.stdout)
# stdlog.setLevel(logging.INFO)
# formatter = logging.Formatter('%(asctime)-12s: %(levelname)-8s\n    %(message)s')
# stdlog.setFormatter(formatter)
# logger.addHandler(stdlog)


class Files(object):
    '''
        open file as file handle for iterator
    '''

    def __init__(self, File):
        self._fos = File

    def __iter__(self):
        if self._fos.lower().endswith(('.gz', '.gzip')):
            fgz = True
            fp = gzip.open(self._fos)
        elif self._fos.lower().endswith('.bz2'):
            fgz = True
            fp = bz2.BZ2File(self._fos)
        else:
            fgz = False
            fp = open(self._fos)

        for index, line in enumerate(fp):
            if index % 50000 == 0 and index != 0:
                logger.info('working on line NO. of {}: {:,}'.format(
                    os.path.basename(self._fos), index
                ))
            try:
                if fgz:
                    line = line.decode()
                # if isinstance(line, bytes):
                # line = line.decode('utf-8')
            except:
                pass
            yield line
        fp.close()


class fpCollapse(object):
    def __init__():
        pass
