#!/usr/bin/env python

import sys

try:
    _, tre_in, name, tre_out = sys.argv
except:
    print('USAGE:\n    python script <tre.nwk> <name.lst> <tre_out.nwk>')
    print('\nNOTES: species rename of nwk phylogentic tree')
    sys.exit(1)


def rename(tre, names, output):
    '''
    phylogentic tree rename
    '''
    with open(tre) as fp1, open(names) as fp2:
        tre = fp1.read()
        names = [i.strip().split('\t')[:2] for i in fp2 if i.strip()]
        names = {i[0]:i[1] for i in names}

    newtre, spid = '', ''
    for i in tre:
        if i in [',', '(', ')', ':', '\n']:
            if spid:
                newtre += names.get(spid, spid)
                spid = ''
            newtre += i
        else:
            spid += i

    with open(output, 'w') as fp:
        fp.write(newtre)


if __name__ == '__main__':
    rename(tre_in, name, tre_out)