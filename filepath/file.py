#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import bz2
import gzip
import time


def LinkFile(source, target, relative=False, backup=True):
    '''
    Create a symbolic link({target}) pointing to {source}.

    Parameters
    ----------
    source:   string, file of absolute path
    basedir:  string
    relative: bool
    backup:   bool, whether backup if target is exists
    '''
    assert os.path.exists(source)
    source = os.path.abspath(source)
    target = os.path.abspath(target)
    if os.path.isdir(target):
        target = '{}{}{}'.format(target, os.sep, os.path.basename(source))
    if os.path.exists(target):
        if os.path.islink(target) and os.path.realpath(target) == source:
            print('file link is exists and right, skip!')
            return
        if backup:
            print('link\'s target is another source file, backup link!')
            FileBackup(target)
        else:
            os.remove(target)
    dirname, basename = os.path.split(target)
    if relative:
        source = RelativePath(source, dirname)
    print('link file "" to "{}"'.format(source, target))
    os.symlink(source, target)


def FileBackup(file):
    '''
    '''
    if os.path.exists(file):
        dirname, basename = os.path.split(file)
        count = 1    # time.strftime("%Y-%m-%d %H:%M:%S")
        timer = time.strftime("%Y%m%d-%H%M")
        backupfile = '{}.{}-{:0>2}.bak'.format(file, timer, count)
        while os.path.exists(backupfile):
            count += 1
            backupfile = '{}.{}-{:0>2}.bak'.format(file, timer, count)
        print("Backing up {0} to {1}".format(file, backupfile))
        os.rename(file, backupfile)
    else:
        print('{0} doesn\'t exist! skip this links ... '.format(file))


def FileName(path):
    sep = os.sep
    fname = path[path.rindex(sep)+1:] if (sep in path) else path
    suffix = {
        'sh', 'shell', 'bash', 'py', 'pl', 'r',
        'gz', 'tar', 'bz', 'bz2', 'xz',
        'fq', 'fa', 'fastq', 'fasta', 'sam', 'bam',
        'bed', 'gtf', 'txt', 'tsv', 'csv',
        'yaml', 'json', 'cfg', 'conf'
    }
    while True:
        index = fname.rindex('.') if ('.' in fname) else -1
        if fname[index+1:].lower() not in suffix:
            return fname
        fname = fname[:index]


def RelativePath(source, basedir):
    '''
    convert abspath to realtive path

    Parameters
    ----------
    source: string
        file or dir of absolute path
    basedir: string
        abspath dir
    '''
    assert os.path.exists(basedir) and os.path.isdir(basedir), 'basedir ...'
    sep = os.sep
    source = os.path.abspath(source).split(sep)
    basedir = os.path.abspath(basedir).split(sep)
    # if not os.path.isdir(basedir):
    #     basedir = os.path.dirname(basedir)
    assert source != basedir, 'source and basedir is same, check!'
    depth = 1
    while True:
        if source[:depth] == basedir[:depth]:
            depth += 1
            continue
        else:
            del source[:(depth-1)]
            del basedir[:(depth-1)]
            break
    basedir = './' if not basedir else '..{}'.format(sep)*len(basedir)
    return os.path.join(basedir, sep.join(source))
