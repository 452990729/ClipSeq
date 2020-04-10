#!/usr/bin/env python2
# coding=utf-8
import sys
import re
import os
import ConfigParser
import argparse

BasePath = os.path.split(os.path.realpath(__file__))[0]
config = ConfigParser.ConfigParser()
config.read(BasePath+'/config.ini')

def HandleConfig(path_in):
    if path_in.startswith('..'):
        return os.path.join(BasePath, path_in)
    else:
        return path_in

#### SOFT
PYTHON = HandleConfig(config.get('SOFTWARE', 'python'))
FASTP = HandleConfig(config.get('SOFTWARE', 'fastp'))
SNAKEMAKE = HandleConfig(config.get('SOFTWARE', 'snakemake'))
STAR = HandleConfig(config.get('SOFTWARE', 'star'))
FASTP = HandleConfig(config.get('SOFTWARE', 'fastp'))
CUTADAPTER = HandleConfig(config.get('SOFTWARE', 'cutadapt'))
SAMTOOLS = HandleConfig(config.get('SOFTWARE', 'samtools'))
BEDTOOLS = HandleConfig(config.get('SOFTWARE', 'bedtools'))
BEDTOBIGBED = HandleConfig(config.get('SOFTWARE', 'bedToBigBed'))
CLIPPER = HandleConfig(config.get('SOFTWARE', 'clipper'))
FASTQSORT = HandleConfig(config.get('SOFTWARE', 'fastq-sort'))
PICARD = HandleConfig(config.get('SOFTWARE', 'picard'))
BEDGRAPHTOBIGWIG = HandleConfig(config.get('SOFTWARE', 'bedGraphToBigWig'))
PERL = HandleConfig(config.get('SOFTWARE', 'perl'))
SNAKEMAKE = HandleConfig(config.get('SOFTWARE', 'snakemake'))
ENV = HandleConfig(config.get('SOFTWARE', 'ENV'))
LD = HandleConfig(config.get('SOFTWARE', 'LD'))

#### SCRIPT
FIXSCORES = HandleConfig(config.get('SCRIPT', 'fix_scores'))
MAKEBIGWIG = HandleConfig(config.get('SCRIPT', 'make_bigwig_files'))
CREATEBIGWIG = HandleConfig(config.get('SCRIPT', 'MakeBigwig'))
COUNTBAM = HandleConfig(config.get('SCRIPT', 'count_aligned_from_sam'))
RMDUPPE = HandleConfig(config.get('SCRIPT', 'barcode_collapse_pe'))
MERGEBAM = HandleConfig(config.get('SCRIPT', 'MergeBam'))
#### DATABASE
HG19 = HandleConfig(config.get('DATABASE', 'hg19'))
HG19SIZES = HandleConfig(config.get('DATABASE', 'hg19sizes'))
REPBASE = HandleConfig(config.get('DATABASE', 'repbase'))

class ReadList(object):
    def __init__(self, line_in):
        list_split = re.split('\t', line_in)
        self.Sample = list_split[0]
        self.Name = list_split[1]
        list_fq = re.split(',', list_split[2])
        self.fq1 = list_fq[0]
        self.fq2 = list_fq[1]
#        self.Genda = list_split[3]

class Snake(object):
    def __init__(self, process):
        self.process = process
        self.input = ''
        self.output = ''
        self.params = ''
        self.log = ''
        self.threads = ''
        self.shell = ''

    def UpdateInput(self, line_in):
        self.input = line_in

    def UpdateOutput(self, line_in):
        self.output = line_in

    def UpdateParams(self, line_in):
        self.params = line_in

    def UpdateLog(self, line_in):
        self.log = line_in

    def UpdateThreads(self, line_in):
        self.threads = line_in

    def UpdateShell(self, line_in):
        self.shell = line_in

    def WriteStr(self, fn):
        fn.write('rule '+self.process+':\n')
        fn.write('\tinput:\n\t\t'+self.input+'\n')
        if self.output:
            fn.write('\toutput:\n\t\t'+self.output+'\n')
        if self.params:
            fn.write('\tparams:\n\t\t'+self.params+'\n')
        if self.log:
            fn.write('\tlog:\n\t\t'+self.log+'\n')
        if self.threads:
            fn.write('\tthreads: '+self.threads+'\n')
        if self.shell:
            fn.write('\tshell:\n\t\t'+self.shell+'\n')
        fn.write('\n')

def MakeSnake(file_in, file_out, platform, cores):
    out = open(file_out, 'w')
    out.write('export PATH={}:$PATH\n'.format(ENV))
    out.write('export LD_LIBRARY_PATH={}:$LD_LIBRARY_PATH\n'.format(LD))
    if platform == 'SGE':
        out.write(SNAKEMAKE+' --cluster "qsub -cwd -V -b y -S /bin/bash -o {log.o} -e {log.e}" -j '+cores+' -s '+file_in+' --printshellcmds --latency-wait 10')
    elif platform == 'local':
        out.write('{} -j {} -s {} --printshellcmds --latency-wait 10'\
                 .format(SNAKEMAKE, cores, file_in))
    out.close()

def main():
    parser = argparse.ArgumentParser(description="WES pipeline")
    parser.add_argument('-c', help='the input fastq list', required=True)
    parser.add_argument('-o', help='the output path', required=True)
    parser.add_argument('-p', help='the platform', choices=['SGE', 'local'], default='local')
    parser.add_argument('-j', help='the job parallel cores', default='2')
    parser.add_argument('-r', help='run now', action='store_true')
    argv=vars(parser.parse_args())
    outpath = argv['o']
    snakefile = open(os.path.join(outpath, 'snakefile.txt'), 'w')
    RawData = os.path.join(outpath, 'RawData')
    list_ob = []
    if not os.path.exists(RawData):
        os.mkdir(RawData)
    else:
        os.system('rm -rf '+RawData)
        os.mkdir(RawData)
    with open(argv['c'], 'r') as f:
        os.chdir(RawData)
        for line in f:
            if not line.startswith('#'):
                ob = ReadList(line.strip())
                list_ob.append(ob)
                os.system('ln -s {} {}'.format(ob.fq1, ob.Sample+'_1.fq.gz'))
                os.system('ln -s {} {}'.format(ob.fq2, ob.Sample+'_2.fq.gz'))
    os.chdir(outpath)
    ### config file
    snakefile.write('Samples = "{}".split()\n'.format(' '.join([i.Sample for i in\
                                                              list_ob])))
    snakefile.write('Groups = "{}".split()\n'.format(' '.join([i.Name for i in\
                                                                list_ob])))


    ###All
    All = Snake('All')
    All.UpdateInput('expand("6.Peak/{group}.final.r2.norm.pos.bw", group=Groups), expand("6.Peak/{group}.final.r2.norm.neg.bw", group=Groups), expand("6.Peak/{group}.final.r2.peaks.fixed.bb", group=Groups)')
    All.WriteStr(snakefile)
    ### QC
    QC = Snake('QC')
    QC.UpdateInput('a = "RawData/{sample}_1.fq.gz", b = "RawData/{sample}_2.fq.gz"')
    QC.UpdateOutput('a = "1.QC/{sample}_1.clean.fq.gz", b = "1.QC/{sample}_2.clean.fq.gz", c = "1.QC/{sample}_QC_report.json", d = "1.QC/{sample}_QC_report.html"')
    QC.UpdateLog('e = "logs/{sample}.qc.e", o = "logs/{sample}.qc.o"')
    QC.UpdateShell(r'"'+FASTP+r' -i {input.a} -o {output.a} -I {input.b} -O {output.b} -A -j {output.c} -h {output.d}"')
    QC.WriteStr(snakefile)
    ###Cutadapter
    Cutadapter = Snake('Cutadapter')
    Cutadapter.UpdateInput('a = "1.QC/{sample}_1.clean.fq.gz", b = "1.QC/{sample}_2.clean.fq.gz"')
    Cutadapter.UpdateOutput('a = "2.Trim/{sample}_1.adapterTrim.fastq.gz", b = "2.Trim/{sample}_2.adapterTrim.fastq.gz"')
    Cutadapter.UpdateLog('e = "logs/Cutadapter.e", o = "logs/Cutadapter.o"')
    Cutadapter.UpdateShell(r'"'+CUTADAPTER+' -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o {output.a} -p {output.b} {input.a} {input.b}"')
    Cutadapter.WriteStr(snakefile)
    ### Cutadapter2
    Cutadapter2 = Snake('Cutadapter2')
    Cutadapter2.UpdateInput('a = "2.Trim/{sample}_1.adapterTrim.fastq.gz", b = "2.Trim/{sample}_2.adapterTrim.fastq.gz"')
    Cutadapter2.UpdateOutput('a = "2.Trim/{sample}_1.adapterTrim_round2.fastq.gz", b = "2.Trim/{sample}_2.adapterTrim_round2.fastq.gz"')
    Cutadapter2.UpdateLog('e = "logs/Cutadapter2.e", o = "logs/Cutadapter2.o"')
    Cutadapter2.UpdateShell(r'"'+CUTADAPTER+' -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 -m 18 -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o {output.a} -p {output.b} {input.a} {input.b}"')
    Cutadapter2.WriteStr(snakefile)
    ### RmRep
    RmRep = Snake('RmRep')
    RmRep.UpdateInput('a = "2.Trim/{sample}_1.adapterTrim_round2.fastq.gz", b = "2.Trim/{sample}_2.adapterTrim_round2.fastq.gz"')
    RmRep.UpdateOutput('a = "3.RmRep/{sample}.rep.bam", b = "3.RmRep/{sample}.rep.bamUnmapped.out.mate1", c = "3.RmRep/{sample}.rep.bamUnmapped.out.mate2"')
    RmRep.UpdateLog('e = "logs/RmRep.e", o = "logs/RmRep.o"')
    RmRep.UpdateThreads('4')
    RmRep.UpdateShell(r'"'+STAR+' --runMode alignReads --runThreadN {threads} --genomeDir '+REPBASE+' --readFilesIn {input.a} {input.b} --outSAMunmapped Within --outFilterMultimapNmax 30 --outFilterMultimapScoreRange 1 --outFileNamePrefix {output.a} --outSAMattributes All --readFilesCommand zcat --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo SM:illumina --alignEndsType EndToEnd > {output.a}"')
    RmRep.WriteStr(snakefile)
    ### StatRep
    StatRep = Snake('StatRep')
    StatRep.UpdateInput('"3.RmRep/{sample}.rep.bam"')
    StatRep.UpdateOutput('"3.RmRep/{sample}.rep.bam.metrics"')
    StatRep.UpdateLog('e = "logs/StatRep.e", o = "logs/StatRep.o"')
    StatRep.UpdateShell(r'"'+SAMTOOLS+' view {input} | '+COUNTBAM+' > {output}"')
    StatRep.WriteStr(snakefile)
    ### FastqSort
    FastqSort = Snake('FastqSort')
    FastqSort.UpdateInput('"3.RmRep/{sample}.rep.bamUnmapped.out.mate1"')
    FastqSort.UpdateOutput('"3.RmRep/{sample}.rep.bamUnmapped.out.sorted.mate1"')
    FastqSort.UpdateLog('e = "logs/FastqSort.e", o = "logs/FastqSort.o"')
    FastqSort.UpdateShell(r'"'+FASTQSORT+' --id {input} > {output}"')
    FastqSort.WriteStr(snakefile)
    ### FastqSort2
    FastqSort2 = Snake('FastqSort2')
    FastqSort2.UpdateInput('"3.RmRep/{sample}.rep.bamUnmapped.out.mate2"')
    FastqSort2.UpdateOutput('"3.RmRep/{sample}.rep.bamUnmapped.out.sorted.mate2"')
    FastqSort2.UpdateLog('e = "logs/FastqSort2.e", o = "logs/FastqSort2.o"')
    FastqSort2.UpdateShell(r'"'+FASTQSORT+' --id {input} > {output}"')
    FastqSort2.WriteStr(snakefile)
    ### Align
    Align = Snake('Align')
    Align.UpdateInput('a = "3.RmRep/{sample}.rep.bamUnmapped.out.sorted.mate1", b = "3.RmRep/{sample}.rep.bamUnmapped.out.sorted.mate2"')
    Align.UpdateOutput('"4.Align/{sample}.rep.align.bam"')
    Align.UpdateLog('e = "logs/Align.e", o = "logs/Align.o"')
    RmRep.UpdateThreads('4')
    Align.UpdateShell(r'"'+STAR+' --runMode alignReads --runThreadN {threads} --genomeDir '+HG19+' --readFilesIn {input.a} {input.b} --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 --outFileNamePrefix {output} --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo SM:illumina --alignEndsType EndToEnd > {output}"')
    Align.WriteStr(snakefile)
    ### RmDup
    RmDup = Snake('RmDup')
    RmDup.UpdateInput('"4.Align/{sample}.rep.align.bam"')
    RmDup.UpdateOutput('a = "4.Align/{sample}.rep.align.rmDup.bam", b = "4.Align/{sample}.rep.align.rmDup.metrics"')
    RmDup.UpdateLog('e = "logs/RmDup.e", o = "logs/RmDup.o"')
    RmDup.UpdateShell(r'"'+RMDUPPE+' --bam {input} --out_file {output.a} --metrics_file {output.b}"')
    RmDup.WriteStr(snakefile)
    ### SortSam
    SortSam = Snake('SortSam')
    SortSam.UpdateInput('"4.Align/{sample}.rep.align.rmDup.bam"')
    SortSam.UpdateOutput('"4.Align/{sample}.rep.align.rmDup.sorted.bam"')
    SortSam.UpdateLog('e = "logs/SortSam.e", o = "logs/SortSam.o"')
    SortSam.UpdateShell(r'"' +PICARD+' SortSam INPUT={input} OUTPUT={output} SORT_ORDER=coordinate"')
    SortSam.WriteStr(snakefile)
    ### Merge
    Merge = Snake('Merge')
    Merge.UpdateInput('expand("4.Align/{sample}.rep.align.rmDup.sorted.bam", sample=Samples)')
    Merge.UpdateOutput('a = "5.Merge", b = expand("5.Merge/{group}.final.bam", group=Groups)')
    Merge.UpdateLog('e = "logs/Merge.e", o = "logs/Merge.o"')
    Merge.UpdateShell(r'"' +MERGEBAM+' -i '+argv['c']+' -o {output.a}"')
    Merge.WriteStr(snakefile)
    ### Index
    Index = Snake('Index')
    Index.UpdateInput('"5.Merge/{group}.final.bam"')
    Index.UpdateOutput('"5.Merge/{group}.final.bam.bai"')
    Index.UpdateLog('e = "logs/Index.e", o = "logs/Index.o"')
    Index.UpdateShell(r'"' +SAMTOOLS+' index {input} {output}"')
    Index.WriteStr(snakefile)
    ### GetRead2
    GetRead2 = Snake('GetRead2')
    GetRead2.UpdateInput('a = "5.Merge/{group}.final.bam", b = "5.Merge/{group}.final.bam.bai"')
    GetRead2.UpdateOutput('"5.Merge/{group}.final.r2.bam"')
    GetRead2.UpdateLog('e = "logs/GetRead2.e", o = "logs/GetRead2.o"')
    GetRead2.UpdateShell(r'"' +SAMTOOLS+' view -hb -f 128 {input.a} > {output}"')
    GetRead2.WriteStr(snakefile)
    ### MakeBigwig
    MakeBigwig = Snake('MakeBigwig')
    MakeBigwig.UpdateInput('"5.Merge/{group}.final.r2.bam"')
    MakeBigwig.UpdateOutput('a = "6.Peak/{group}.final.r2.norm.pos.bw", b = "6.Peak/{group}.final.r2.norm.neg.bw"')
    MakeBigwig.UpdateLog('e = "logs/MakeBigwig.e", o = "logs/MakeBigwig.o"')
    MakeBigwig.UpdateShell(r'"' +CREATEBIGWIG+' --bam {input} --genome '+HG19SIZES+' --bw_pos {output.a} --bw_neg {output.b}"')
    MakeBigwig.WriteStr(snakefile)
    ### CallPeak
    CallPeak = Snake('CallPeak')
    CallPeak.UpdateInput('"5.Merge/{group}.final.r2.bam"')
    CallPeak.UpdateOutput('"6.Peak/{group}.final.r2.peaks.bed"')
    CallPeak.UpdateLog('e = "logs/CallPeak.e", o = "logs/CallPeak.o"')
    CallPeak.UpdateThreads('4')
    CallPeak.UpdateShell(r'"' +CLIPPER+' -b {input} -s hg19 -o {output} --save-pickle --processors {threads}"')
    CallPeak.WriteStr(snakefile)
    ### FixPeak
    FixPeak = Snake('FixPeak')
    FixPeak.UpdateInput('"6.Peak/{group}.final.r2.peaks.bed"')
    FixPeak.UpdateOutput('"6.Peak/{group}.final.r2.peaks.fixed.bed"')
    FixPeak.UpdateLog('e = "logs/FixPeak.e", o = "logs/FixPeak.o"')
    FixPeak.UpdateShell(r'"' +FIXSCORES+' --bed {input} --out_file {output}"')
    FixPeak.WriteStr(snakefile)
    ### BedToBigBed
    BedToBigBed = Snake('BedToBigBed')
    BedToBigBed.UpdateInput('"6.Peak/{group}.final.r2.peaks.fixed.bed"')
    BedToBigBed.UpdateOutput('"6.Peak/{group}.final.r2.peaks.fixed.bb"')
    BedToBigBed.UpdateLog('e = "logs/BedToBigBed.e", o = "logs/BedToBigBed.o"')
    BedToBigBed.UpdateShell(r'"'+BEDTOBIGBED+' {input} '+HG19SIZES+' {output} -type=bed6+4"')
    BedToBigBed.WriteStr(snakefile)

    ######RUN
    OutShell = os.path.join(outpath, 'work.sh')
    MakeSnake(os.path.join(outpath, 'snakefile.txt'), OutShell, argv['p'], argv['j'])
    if argv['r']:
        os.system('nohup sh {}&'.format(OutShell))

if __name__ == '__main__':
    main()
