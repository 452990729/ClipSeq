#!/usr/bin/env python2


import sys
import os
import re
import argparse
import ConfigParser

BasePath = os.path.split(os.path.realpath(__file__))[0]
config = ConfigParser.ConfigParser()
config.read(BasePath+'/../Bin/config.ini')

def HandleConfig(path_in):
    if path_in.startswith('..'):
        return os.path.join(BasePath+'/../Bin/', path_in)
    else:
        return path_in

SAMTOOLS = HandleConfig(config.get('SOFTWARE', 'samtools'))

def ReadTable(file_in):
    dict_g = {}
    with open(file_in, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                list_split = re.split('\t', line.strip())
                if list_split[1] not in dict_g:
                    dict_g[list_split[1]] = [list_split[0]]
                else:
                    dict_g[list_split[1]] += [list_split[0]]
    return dict_g

def MergeFile(dict_g, outpath):
    for g in dict_g:
        list_tmp = ['4.Align/{}.rep.align.rmDup.sorted.bam'.format(i) for i in dict_g[g]]
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        out = '{}/{}.final.bam'.format(outpath, g)
        os.system('{} merge {} {}'.format(SAMTOOLS, out, ' '.join(list_tmp)))

def main():
    parser = argparse.ArgumentParser(description="Merge Bam By Group")
    parser.add_argument('-i', help='the input fastq list', required=True)
    parser.add_argument('-o', help='output path', required=True)
    argv=vars(parser.parse_args())
    dict_g = ReadTable(argv['i'])
    MergeFile(dict_g, argv['o'])


if __name__ == '__main__':
    main()

