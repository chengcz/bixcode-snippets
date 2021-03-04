#!/usr/bin/env python
# -*- coding:utf-8 -*-


def intervals(start, end):
    '''
    parameter
    ---------
    start:    int/list of int
    end:      int/list of int
    '''
    if isinstance(start, list) and isinstance(end, list):
        return Interval(*[AtomicInterval(x, y) for x, y in zip(start, end)])
    elif isinstance(start, int) and isinstance(end, int):
        return Interval(AtomicInterval(start, end))
    else:
        raise Exception('wrong input')


class AtomicInterval(object):
    '''
    '''
    __slots__ = ('_lower', '_upper')
    def __init__(self, start, end):
        self._lower = start
        self._upper = end
        if self.is_empty():
            self._lower = float('inf')
            self._upper = float('-inf')

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f'AtomicInterval: [{self._lower}, {self._upper})'

    @property
    def lower(self):
        return self._lower

    @property
    def upper(self):
        return self._upper

    def is_empty(self):
        return (self._lower >= self._upper)

    def overlaps(self, other, adjacent=False):
        if not isinstance(other, AtomicInterval):
            raise TypeError('Only AtomicInterval instances are supported.')
        if other.is_empty() or self.is_empty():
            return False
        first, second = (self, other) if (self.lower <= other.lower) else (other, self)
        if adjacent:
            return first.upper >= second.lower
        return first.upper > second.lower

    def __len__(self):
        return self.upper - self.lower

    def __eq__(self, other):
        if isinstance(other, AtomicInterval):
            return (self.lower == other.lower) and (self.upper == other.upper)
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.upper <= other.lower
        elif isinstance(other, int):
            return self.upper <= other
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.upper <= other.upper
        elif isinstance(other, int):
            return self.upper <= other
        else:
            return NotImplemented

    def __gt__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.lower >= other.upper
        elif isinstance(other, int):
            return self.lower >= other
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.lower >= other.lower
        elif isinstance(other, int):
            return self.lower >= other
        else:
            return NotImplemented

    def __and__(self, other):
        if isinstance(other, AtomicInterval):
            if self.overlaps(other):
                lower = max(self.lower, other.lower)
                upper = min(self.upper, other.upper)
                return AtomicInterval(lower, upper)
            else:
                return AtomicInterval(float('inf'), float('-inf'))
        else:
            return NotImplemented

    def __or__(self, other):
        if isinstance(other, AtomicInterval):
            if self.overlaps(other, adjacent=True):
                lower = min(self.lower, other.lower)
                upper = max(self.upper, other.upper)
                return AtomicInterval(lower, upper)
            else:
                return Interval(self, other)
        else:
            return NotImplemented

    def __contains__(self, item):
        if isinstance(item, AtomicInterval):
            left  = self.lower <= item.lower
            right = item.upper <= self.upper
            return left and right
        elif isinstance(item, Interval):
            for interval in item:
                if interval not in self:
                    return False
            return True
        elif isinstance(item, int):
            return self.lower <= item < self.upper
        else:
            return NotImplemented

    def __invert__(self):
        if self.is_empty():
            return AtomicInterval(float('-inf'), float('inf'))
        return Interval(
            AtomicInterval(float('-inf'), self.lower),
            AtomicInterval(self.upper, float('inf'))
        )

    def __sub__(self, other):
        if isinstance(other, AtomicInterval):
            return self & ~other
        else:
            return NotImplemented


class Interval(object):
    '''
    '''
    __slots__ = ('_intervals')
    def __init__(self, *intervals):
        self._intervals = list()

        for x in intervals:
            if isinstance(x, Interval):
                self._intervals.extend(x)
            elif isinstance(x, AtomicInterval):
                # if not x.is_empty():
                    self._intervals.append(x)
            else:
                raise TypeError('Parameters must be Interval or AtomicInterval instances')
        self.__deldup()

    def __deldup(self):
        self._intervals = [x for x in self._intervals if not x.is_empty()]
        if len(self._intervals) == 0:
            self._intervals.append(AtomicInterval(float('inf'), float('-inf')))
        else:
            self._intervals.sort(key=lambda x: (x.lower))

            x = 0
            while x < len(self._intervals) - 1:
                first = self._intervals[x]
                second = self._intervals[x + 1]
                if first.overlaps(second, adjacent=True):
                    interval = first | second
                    self._intervals.pop(x)  # pop first
                    self._intervals.pop(x)  # pop second
                    self._intervals.insert(x, interval)
                else:
                    x = x + 1

    def __getitem__(self, item):
        return self._intervals[item]

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        regions = [f'[{x.lower}, {x.upper})' for x in self]
        return f'Interval: {" / ".join(regions)}'

    @property
    def lower(self):
        return self._intervals[0].lower

    @property
    def upper(self):
        return self._intervals[-1].upper

    def overlaps(self, other):
        if isinstance(other, Interval):
            for x in other._intervals:
                if self.overlaps(x):
                    return True
            return False
        elif isinstance(other, AtomicInterval):
            for x in self._intervals:
                if x.overlaps(other):
                    return True
            return False
        elif isinstance(other, int):
            return other in self
        else:
            return NotImplemented

    def to_atomic(self):
        return Interval(AtomicInterval(self.lower, self.upper))

    def is_atomic(self):
        return len(self._intervals) == 1

    def is_empty(self):
        return self.is_atomic() and self._intervals[0].is_empty()

    def append(self, other):
        if isinstance(other, AtomicInterval):
            self._intervals.append(other)
        elif isinstance(other, Interval):
            self._intervals.extend(other._intervals)
        else:
           raise TypeError('Parameters must be Interval/AtomicInterval instances')
        self.__deldup()

    @property
    def length(self):
        return sum([len(x) for x in self])

    def __len__(self):
        return len(self._intervals)

    def __iter__(self):
        return iter(self._intervals)

    def __eq__(self, other):
        if isinstance(other, Interval):
            return self._intervals == other._intervals
        elif isinstance(other, AtomicInterval):
            return Interval(other) == self
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.upper <= other.lower
        elif isinstance(other, int):
            return self.upper <= other
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.upper <= other.upper
        elif isinstance(other, int):
            return self.upper <= other
        else:
            return NotImplemented

    def __gt__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.lower >= other.upper
        elif isinstance(other, int):
            return self.lower >= other
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.lower >= other.lower
        elif isinstance(other, int):
            return self.lower >= other
        else:
            return NotImplemented

    def __and__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            if isinstance(other, AtomicInterval):
                lst = [other]
            else:
                lst = list(other._intervals)
            new_intervals = []
            for x in self._intervals:
                for y in lst:
                    new_intervals.append(x & y)
            return Interval(*new_intervals)
        else:
            return NotImplemented

    def __rand__(self, other):
        return self & other

    def __or__(self, other):
        if isinstance(other, AtomicInterval):
            return self | Interval(other)
        elif isinstance(other, Interval):
            return Interval(*(self._intervals + other._intervals))
        else:
            return NotImplemented

    def __ror__(self, other):
        return self | other

    def __sub__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self & ~other
        else:
            return NotImplemented

    def __rsub__(self, other):
        return other & ~self

    def __invert__(self):
        complements = [~x for x in self._intervals]
        intersection = complements[0]
        for interval in complements:
            intersection = intersection & interval
        return intersection

    def __contains__(self, item):
        if isinstance(item, Interval):
            for x in item._intervals:
                if x not in self:
                    return False
            return True
        elif isinstance(item, AtomicInterval):
            for x in self._intervals:
                if item in x:
                    return True
            return False
        elif isinstance(item, int):
            for x in self._intervals:
                if item in x:
                    return True
            return False
        else:
            return NotImplemented


def genomic_intervals(chrom, start, end):
    if isinstance(start, int) and isinstance(end, int):
        region = [[start, end],]
    elif isinstance(start, list) and isinstance(end, list):
        region = [
            [x, y] for x, y in zip(start, end)
            if isinstance(x, int) and isinstance(y, int)
        ]
    else:
        raise TypeError('Only list/intger instances ars supported.')
    return GenomicInterval(chrom, Interval(*[AtomicInterval(x, y) for x, y in region]))


class GenomicInterval(object):
    '''
    parameter
    ---------
    chrom:     string, [required]
    start:     int/list, [required]
    end:       int/list, [required]
    name:      string
    strand:    string, (-, +, None)
    '''
    __slots__ = ('_chrom', '_intervals', '_attr')
    def __init__(self, chrom, start, end, name=None, strand=None, **attr):
        self._chrom = chrom
        self._intervals = Interval(*interval)
        self._attr = attr

    @staticmethod
    def __input_to_interval(start, end):
        region = []
        if isinstance(start, int) and isinstance(end, int):
            region.append(AtomicInterval(start, end))
        elif isinstance(start, list) and isinstance(end, list):
            if len(start) != len(end):
                raise Exception('length of start and end is unequal in Input interval')
            for x, y in zip(start, end):
                if not (isinstance(x, int) and isinstance(y, int)):
                    raise TypeError('Input interval data type is not int')
                region.append(AtomicInterval(x, y))
        else:
            raise TypeError('Unsupported input data type!')
        return Interval(region)

    def __str__(self):
        if self.strand:
            return f'GenomicInterval: {self.chrom}:{self.start}-{self.end}:{self.strand}'
        else:
            return f'GenomicInterval: {self.chrom}:{self.start}-{self.end}'

    def __repr__(self):
        return self.__str__()

    @property
    def chrom(self):
        return self._chrom

    @property
    def start(self):
        return self._chrom

    @property
    def end(self):
        return self._region

    @property
    def strand(self):
        return self._chrom

    @property
    def name(self):
        return self._region

    @property
    def interval(self):
        return self._intervals

    def append_interval(self, start, end):
        region = AtomicInterval(start, end)
        self._region = self._region | region

    def mask_interval(self, start, end):
        region = AtomicInterval(start, end)
        self._region = self._region - region

    def overlaps(self, other):
        if isinstance(other, GenomicInterval):
            if self._chrom == other.chrom:
                if self._region.overlaps(self.region):
                    return True
            return False
        else:
            return NotImplemented

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        regions = [f'{x.lower} - {x.upper}' for x in self._region]
        return f'GenomicInterval: {self._chrom}: {" / ".join(regions)}'

    def __len__(self):
        return sum([len(x) for x in self._region])

    def __eq__(self, other):
        if isinstance(other, GenomicInterval):
            if (self._chrom == other.chrom) and (self._region == other.region):
                return True
            return False
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if isinstance(other, GenomicInterval):
            if (self._chrom == other.chrom):
                return self._region < other.region
            return False
        elif isinstance(other, int):
            return self._region.upper < other
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, GenomicInterval):
            if (self._chrom == other.chrom):
                return self._region <= other.region
            return False
        elif isinstance(other, int):
            return self._region.upper <= other
        else:
            return NotImplemented

    def __gt__(self, other):
        if isinstance(other, GenomicInterval):
            if (self._chrom == other.chrom):
                return self._region > other.region
            return False
        elif isinstance(other, int):
            return self._region.lower > other
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, GenomicInterval):
            if (self._chrom == other.chrom):
                return self._region >= other.region
            return False
        elif isinstance(other, int):
            return self._region.lower >= other
        else:
            return NotImplemented

    def __and__(self, other):
        if isinstance(other, GenomicInterval):
            if self._chrom != other.chrom:
                return NotImplemented
            new_intervals = []
            for interval in self.region:
                for o_interval in other.region:
                    new_intervals.append(interval & o_interval)
            return GenomicInterval(self.chrom, Interval(*new_intervals))
        else:
            return NotImplemented

    def __rand__(self, other):
        return self & other

    def __or__(self, other):
        if isinstance(other, GenomicInterval):
            if self._chrom != other.chrom:
                return NotImplemented
            region = Interval(*(self.region._intervals + other.region._intervals))
            return GenomicInterval(self.chrom, region)
        else:
            return NotImplemented

    def __ror__(self, other):
        return self | other

    def __sub__(self, other):
        if isinstance(other, GenomicInterval):
            if self._chrom != other.chrom:
                return NotImplemented
            return self & ~other
        else:
            return NotImplemented

    def __rsub__(self, other):
        return other & ~self

    def __invert__(self):
        intervals = [~x for x in self.region._intervals]
        intersection = intervals[0]
        for interval in intervals:
            intersection = intersection & interval
        return GenomicInterval(self.chrom, Interval(*intersection))

    def __contains__(self, other):
        if isinstance(other, GenomicInterval):
            if self._chrom != other.chrom:
                return False
            for o_interval in other.region._intervals:
                if o_interval not in self.region:
                    return False
            return True
        elif isinstance(other, int):
            for interval in self.region._intervals:
                if other in interval:
                    return True
            return False
        else:
            return NotImplemented

