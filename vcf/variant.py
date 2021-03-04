

class Header(object):
    def __init__(self):
        pass


class genotype(object):
    def __init__(self):
        pass


class Variant(object):
    '''
    parameter
    ---------
    chrom:    string
    pos:      integer, 1-based location system
    ref:      string
    alt:      string
    attr:     dict
    '''
    __slots__ = ('_chrom', '_pos', '_ref', '_alt', '_attr', '_varid')
    def __init__(self, chrom: str, pos: int, ref: str, alt: str, attri: dict = None):
        self._chrom = self.__chrom_modifer(chrom)
        self._pos = pos
        self._ref = ref
        self._alt = alt
        self._varid = self.__generate_varid()
        attri = attri if isinstance(attri, dict) else {}
        self._attr = attri

    @staticmethod
    def __chrom_modifer(chrom):
        if chrom.lower().startswith('chr') and ('_' not in chrom):
            chrom = chrom[3:]
        if chrom.lower() == 'm':
            chrom == 'MT'
        return chrom.upper()

    def __generate_varid(self):
        varid = f'{self._chrom}_{self._pos}_{self._ref}_{self._alt}'
        varid = hashlib.md5(varid.encode('utf-8')).hexdigest()
        return varid

    def __str__(self):
        return f'<Variant object: {self._chrom}:{self._pos}{self._ref}>{self._alt} ... >'

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if not (self._pos, self._ref, self._alt) == (other.pos, other.ref, other.alt):
            return False
        if not (self._chrom == other.chrom):
            return False
        return True

    @property
    def varid(self):
        return self._varid

    @property
    def chrom(self):
        return self._chrom

    @property
    def pos(self):
        return self._pos

    @property
    def ref(self):
        return self._ref

    @property
    def alt(self):
        return self._alt

    @property
    def id(self):
        return self._attr.get('ID', None)

    @property
    def start(self):
        return self._pos - 1

    @property
    def end(self):
        return self._pos + len(self._ref) - 1

    @property
    def info(self, parse=True):
        if parse:
            return self.__parse_INFO()
        else:
            return self._attr.get('INFO', '.')

    def __parse_INFO(self):
        infor = self._attr.get('INFO', None)
        if infor:
            infor = [x.split('=', 1) for x in infor.split(';')]
            info = {}
            for x in infor:
                if len(x) == 2:
                    info[x[0]] = x[1]
                elif len(x) == 1:
                    info[x] = True
            # info = {x[0]: x[1] for x in infor if len(x) == 2}
            return info
        else:
            return {}

    def return_vcf(self, fill_sample=False):
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
        if fill_sample:
            return (f'{self._chrom}\t{self._pos}\t.\t{self._ref}\t{self._alt}'
                    '\t.\t.\t.\tGT:AD:DP:GQ:PL\t1/1:.:.:.:.\n')
        else:
            return f'{self._chrom}\t{self._pos}\t.\t{self._ref}\t{self._alt}\t.\t.\t.\n'
