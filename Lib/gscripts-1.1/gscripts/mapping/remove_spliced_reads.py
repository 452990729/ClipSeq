#!/usr/bin/env python


__author__ = 'olga'
from glob import glob
from gscripts.qtools._Submitter import Submitter
import sys
import argparse


class CommandLine(object):
    """
    Check out the argparse documentation for more awesome things you can do!
    https://docs.python.org/2/library/argparse.html
    """

    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Removed reads that are spliced from a bam file')
        parser.add_argument('job_name', required=False,
                            type=str, action='store', default='remove_spliced',
                            help='Sample info file with bam files and sample '
                                 'IDs')
        parser.add_argument('--do-not-submit', required=False,
                            action='store_true',
                            help='Flag to not actually submit the job but '
                                 'just write the sh file (for testing)')
        parser.add_argument('-o', '--out-sh', required=False, type=str,
                            action='store',
                            help='Name of the sh file to submit to the job '
                                 'scheduler')
        parser.add_argument('--queue-type', required=False, type=str,
                            action='store', default='PBS',
                            help='Type of the queue to submit to. For testing '
                                 'purposes on non-server machines (e.g. '
                                 'laptops)')
        parser.add_argument('-d', '--directory', required=False,
                            action='store', default='./',
                            help='Directory where the bam files are. Default '
                                 'is the current directory.')

        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))

    def do_usage_and_die(self, str):
        '''
        If a critical error is encountered, where it is suspected that the
        program is not being called with consistent parameters or data, this
        method will write out an error string (str), then terminate execution
        of the program.
        '''
        import sys

        print >> sys.stderr, str
        self.parser.print_usage()
        return 2


# Class: Usage
class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual
    exit when raised
    '''

    def __init__(self, msg):
        self.msg = msg


class RemoveSpliced(object):
    def __init__(self, job_name, out_sh=None,
                 directory='./', submit=True):
        cmd_list = []
        for file in glob('{}/*bam'.format(directory.rstrip('/'))):
            cmd_list.append(
                "samtools view -h -F 4 {0} | awk '$6 !~ /N/ || $1 ~ /@/' "
                "| "
                "samtools view -bS - > {0}.unspliced.bam".format(file))

        sub = Submitter(queue_type='PBS', sh_filename=out_sh,
                        commands=cmd_list,
                        job_name=job_name, nodes=1, ppn=16, queue='home',
                        array=True,
                        walltime='0:30:00',
                        max_running=10)
        sub.job(submit=submit)


if __name__ == '__main__':
    cl = CommandLine()

    job_name = cl.args['job_name']
    out_sh = job_name + '.sh' if cl.args['out_sh'] is None \
        else cl.args['out_sh']
    submit = not cl.args['do_not_submit']
    directory = cl.args['directory'].rstrip('/')

    RemoveSpliced(job_name, out_sh, directory, submit)
