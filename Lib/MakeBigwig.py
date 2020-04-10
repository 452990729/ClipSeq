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

MAKEBIGWIG = HandleConfig(config.get('SCRIPT', 'make_bigwig_files'))
BEDGRAPHTOBIGWIG = HandleConfig(config.get('SOFTWARE', 'bedGraphToBigWig'))

def MakeFile(bam, genome, bw_pos, bw_neg):
    shell = '{} --bam {} --genome {} --bw_pos {} --bw_neg {}'.format(MAKEBIGWIG, bam, genome, bw_pos, bw_neg)
    os.system(shell)

def SortFile(infile, outfile):
    shell = 'sort -k1,1 -k2,2n {} > {}'.format(infile, outfile)
    os.system(shell)

def MakeBigWig(infile, outfile, genome):
    shell = '{} {} {} {}'.format(BEDGRAPHTOBIGWIG, infile, genome, outfile)
    os.system(shell)

def Process(bam, genome, bw_pos, bw_neg):
    outpath = os.path.dirname(bw_pos)
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    MakeFile(bam, genome, bw_pos, bw_neg)
    lb = '.'.join(re.split('\.', bam)[:-1])
    SortFile(lb+'.norm.pos.bg', lb+'.norm.pos.sort.bg')
    SortFile(lb+'.norm.neg.t.bg', lb+'.norm.neg.t.sort.bg')
    MakeBigWig(lb+'.norm.pos.sort.bg', bw_pos, genome)
    MakeBigWig(lb+'.norm.neg.t.sort.bg', bw_neg, genome)

def main():
    parser = argparse.ArgumentParser(description="Merge Bam By Group")
    parser.add_argument('--bam', help='the input bam', required=True)
    parser.add_argument('--genome', help='the genome size file', required=True)
    parser.add_argument('--bw_pos', help='output pos', required=True)
    parser.add_argument('--bw_neg', help='output neg', required=True)
    argv=vars(parser.parse_args())
    Process(argv['bam'], argv['genome'], argv['bw_pos'], argv['bw_neg'])


if __name__ == '__main__':
    main()

