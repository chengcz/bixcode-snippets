#!/usr/bin/env python


class Feature(object):
    '''
    '''
    def __init__(self):
        pass

    @property
    def chrom(self):
        return self._chrom

    @property
    def source(self):
        return self._source

    @property
    def feature(self):
        return self._feature

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def strand(self):
        return self._strand

    @property
    def score(self):
        return self._score

    @property
    def frame(self):
        return self._frame

    @property
    def attrs(self):
        return self._attr


class Transcript(Feature):
    '''
    '''
    def __init__(self, *features):
        pass


class Gene(Feature):
    '''
    '''
    def __init__(self, *transcripts):
        pass
