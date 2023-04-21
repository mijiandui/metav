#!/usr/bin/env python

# -------------------------------------------------------------------------
# VIRSORT SOP
# Copyright (C) 2022 - 2023 Chouxian Ma <biomath_2014@163.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# -------------------------------------------------------------------------


from __future__ import print_function

import getopt
import json
import logging
import multiprocessing
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import math
from common import *

logger = logging.getLogger(__name__)

_usage_message = '''
contact: Chouxian Ma <biomath_2014@163.com>

Usage:
  virsort2_sop [options] {{-i <contig file>}} [-o <out_dir>]

  Input options that can be specified for multiple times (supporting plain text and gz/bz2 extensions)
    -i                       <file>         contig fasta file
    -d/--db-dir              <dir>          checkv database directory

Optional Arguments:
    --include-groups         <str>          classifier of viral groups to include(comma seperated and no space in between);
                                            available options are dsDNAphage,NCLDV,RNA,ssDNA,ssDNA, [dsDNAphage,ssDNA]
    --min-score              <float>        minimal score to be identified as viral, [0.5]
    --min-length             <int>          minimal length required; all seqs shorter than this will be removed [1500]
    -t/--num-cpu-threads     <int>          number of CPU threads [# of logical processors]
    -n                       <int>          number of split [10]

  Output options:
    -o/--out-dir             <string>       output directory [./virsorter2]
    -p/--out-prefix          <string>       output prefix, eg. sample name
    -f/--force                              overwritten exist directory
    --keep-tmp-files                        keep all temporary files
    --tmp-dir                <string>       set temp directory

Other Arguments:
    --no-split                              do not split raw contigs FASTA into multiple parts & run single virsorter instance
    --continue                              continue a METAV run from its last available check point.
                                            please set the output directory correctly when using this option.
    -h/--help                               print the usage message
    -v/--version                            print version
'''


class SoftwareInfo(object):
    script_path = os.path.dirname(os.path.realpath(__file__))
    virsorter2 = 'virsorter'
    checkv = 'checkv'
    parallel_virsort = os.path.join(script_path, 'pvirsort')
    viral_screen_script = os.path.join(script_path, 'viral_screen')

    @property
    def virsorter_sop_version(self):
        return 'VIRSORTER2_SOP v1.0.1'
    
    @property
    def usage_message(self):
        return _usage_message


class Options:
    def __init__(self):
        self.out_dir = ''
        self.temp_dir = ''
        self.continue_mode = False
        self.force_overwrite = False
        self.do_not_split = False
        self.min_length = 1500
        self.num_cpu_threads = 0
        self.keep_tmp_files = False
        self.out_prefix = ''
        self.input_contig_file = ''
        self.db_dir = ''
        self.viral_groups = 'dsDNAphage,ssDNA'
        self.min_score = 0.5
        self.nparts = 10
        self.verbose = False

    @property
    def log_file_name(self):
        if self.out_prefix == '':
            return os.path.join(self.out_dir, 'log')
        else:
            return os.path.join(self.out_dir, self.out_prefix + '.log')

    @property
    def option_file_name(self):
        return os.path.join(self.out_dir, 'options.json')

    @property
    def virsorter_pass1_dir(self):
        return os.path.join(self.out_dir, 'vs2-pass1')

    @property
    def checkv_dir(self):
        return os.path.join(self.out_dir, 'checkv')

    @property
    def checkv_proviruses_contig(self):
        return os.path.join(self.checkv_dir, 'proviruses.fna')

    @property
    def checkv_viruses_contig(self):
        return os.path.join(self.checkv_dir, 'viruses.fna')

    @property
    def checkv_final_contig(self):
        return os.path.join(self.checkv_dir, 'combined.fna')
    
    @property
    def checkv_contamination_table(self):
        return os.path.join(self.checkv_dir, 'contamination.tsv')

    @property
    def virsort2_pass1_contig_file(self):
       return os.path.join(self.virsorter_pass1_dir, 'final-viral-combined.fa')

    @property
    def virsort2_pass1_score_file(self):
       return os.path.join(self.virsorter_pass1_dir, 'final-viral-score.tsv')
     
    @property
    def virsorter2(self):
        return software_info.virsorter2

    @property
    def parallel_virsort(self):
        return software_info.parallel_virsort
    
    @property
    def checkv(self):
        return software_info.checkv

    @property
    def viral_screen_script(self):
        return software_info.viral_screen_script

    def dump(self):
        with open(self.option_file_name, 'w') as f:
            json.dump(self.__dict__, f)

    def load_for_continue(self):
        with open(self.option_file_name, 'r') as f:
            self.__dict__.update(json.load(f))


software_info = SoftwareInfo()
opt = Options()
check_point = Checkpoint(logger)


def run_sub_command(cmd, msg, verbose=False):
    if opt.verbose:
        verbose = opt.verbose

    logger.info(msg)
    cmd_str = ' '.join(cmd)
    logger.info('command %s' % cmd_str)

    start_time = time.time()
    
    try:
        p = subprocess.run(cmd_str, shell=True, check=True, capture_output=True)

        logger.info(p.stderr.decode('utf-8'))

        if p.returncode != 0:
            logger.error('Error occurs, please refer to %s for detail' % opt.log_file_name)
            logger.error('Command: %s; Exit code %d' % (cmd_str, int(ret_code)))
            exit(ret_code)
    except KeyboardInterrupt:
        exit(signal.SIGINT)

    logger.info('CMD DONE. Time elapsed: %f seconds ' % (time.time() - start_time))

def parse_option(argv):
    try:
        opts, args = getopt.getopt(argv, 'hvo:t:i:p:n:d:f',
                                   ['help',
                                    'out-dir=',
                                    'db-dir=',
                                    'num-cpu-threads=',
                                    'keep-tmp-files',
                                    'tmp-dir=',
                                    'include-groups=',
                                    'min-score=',
                                    'min-length=',
                                    'continue',
                                    'version',
                                    'verbose',
                                    'no-split',
                                    'out-prefix=',
                                    'force'])
    except getopt.error as msg:
        raise Usage(software_info.virsorter_sop_version + '\n' + str(msg))
    if len(opts) == 0:
        raise Usage(software_info.virsorter_sop_version + '\n' + software_info.usage_message)

    for option, value in opts:
        if option in ('-h', '--help'):
            print(software_info.virsorter_sop_version + '\n' + software_info.usage_message)
            exit(0)
        elif option in ('-o', '--out-dir'):
            opt.out_dir = abspath(value)
        elif option == '--min-length':
            opt.min_length = int(value)
        elif option in ('-t', '--num-cpu-threads'):
            opt.num_cpu_threads = int(value)
        elif option == '--keep-tmp-files':
            opt.keep_tmp_files = True
        elif option in ('-v', '--version'):
            print(software_info.virsorter_sop_version)
            exit(0)
        elif option == '--verbose':
            opt.verbose = True
        elif option == '--continue':
            opt.continue_mode = True
        elif option == '--no-split':
            opt.do_not_split = True
        elif option in ('-p', '--out-prefix'):
            opt.out_prefix = value
        elif option == '--tmp-dir':
            opt.temp_dir = abspath(value)
        elif option == '-i':
            opt.input_contig_file = value
        elif option == '-n':
            opt.nparts = int(value)
        elif option in ('-f', '--force'):
            opt.force_overwrite = True
        elif option in ('-d', '--db-dir'):
            opt.db_dir = value
        elif option == '--include-groups':
            opt.viral_groups = value
        elif option == '--min-score':
            opt.min_score = float(value)
        else:
            raise Usage('Invalid option {0}'.format(option))


def setup_output_dir():
    if not opt.out_dir:
        opt.out_dir = abspath('./virsorter_checkv_out')

    check_point.set_file(os.path.join(opt.out_dir, 'checkpoints.txt'))

    if opt.continue_mode and not os.path.exists(opt.option_file_name):
        print('Cannot find {0}, switching to normal mode'.format(opt.option_file_name), file=sys.stderr)
        opt.continue_mode = False

    if opt.continue_mode:
        print('Continue mode activated. Ignore all options except for -o/--out-dir.', file=sys.stderr)
        opt.load_for_continue()
        check_point.load_for_continue()
    else:
        if not opt.force_overwrite and os.path.exists(opt.out_dir):
            raise Usage(
                'Output directory ' + opt.out_dir +
                ' already exists, please change the parameter -o to another value to avoid overwriting.')

        if opt.temp_dir == '':
            opt.temp_dir = os.path.join(opt.out_dir, 'tmp')
        else:
            opt.temp_dir = tempfile.mkdtemp(dir=opt.temp_dir, prefix='virsorter2_tmp_')

    mkdir_if_not_exists(opt.out_dir)
    mkdir_if_not_exists(opt.temp_dir)


def setup_logger():
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', '%Y-%m-%d %H:%M:%S')

    file_handler = logging.FileHandler(opt.log_file_name, 'a')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)
    logger.addHandler(console)

    logger.info(software_info.virsorter_sop_version)


def check_and_correct_option():
    # set mode

    if opt.min_length <= 0:
        raise Usage('--min-length must be greater than 0.')
    if opt.min_score < 0 or opt.min_score > 1:
        raise Usage('--min-score must be in range [0, 1]')
    
    if opt.num_cpu_threads > multiprocessing.cpu_count():
        logger.warning('Maximum number of available CPU thread is %d.' % multiprocessing.cpu_count())
        logger.warning('Number of thread is reset to the %d.' % multiprocessing.cpu_count())
        opt.num_cpu_threads = multiprocessing.cpu_count()
    if opt.num_cpu_threads == 0:
        opt.num_cpu_threads = multiprocessing.cpu_count()

    if opt.db_dir == '':
        raise Usage('Please specify db-dir for checkv')
    if not os.path.exists(opt.db_dir):
        raise Usage('db dir dose not exist: ' + opt.db_dir)
    if opt.out_prefix == '':
        opt.out_prefix = 'sample'

def check_input():
    if not os.path.exists(opt.input_contig_file):
        raise Usage('Cannot find file ' + opt.input_contig_file)


@check_point
def virsort_pass1():
    if not opt.continue_mode and os.path.exists(opt.virsorter_pass1_dir):
        shutil.rmtree(opt.virsorter_pass1_dir)

    virsorter_cmd = [opt.virsorter2, 
             'run',
             '--keep-original-seq',
             '-i', opt.input_contig_file, 
             '-w', opt.virsorter_pass1_dir,
             '--min-length', str(opt.min_length),
             '--min-score', str(opt.min_score),
             '--include-groups', opt.viral_groups,
             '--tmpdir', opt.temp_dir,
             '-j', str(opt.num_cpu_threads),
             'all']

    run_sub_command(virsorter_cmd, 'first pass of virsorter2 ')

@check_point
def parallel_virsort_pass1():
    if not opt.continue_mode and os.path.exists(opt.virsorter_pass1_dir):
        shutil.rmtree(opt.virsorter_pass1_dir)

    virsorter_cmd = [opt.parallel_virsort, 
             '-i', opt.input_contig_file, 
             '-o', opt.virsorter_pass1_dir,
             '-l', str(opt.min_length),
             '-s', str(opt.min_score),
             '-g', opt.viral_groups,
             '-n', str(opt.nparts),
             '-t', str(opt.num_cpu_threads),
             'all']

    run_sub_command(virsorter_cmd, 'run virsorter2 in parallel ')


@check_point
def checkv():

    if not opt.continue_mode and os.path.exists(opt.checkv_dir):
        shutil.rmtree(opt.checkv_dir)

    if not os.path.exists(opt.virsort2_pass1_contig_file):
        logger.error('can not find first pass contig file: ' + opt.virsort2_pass1_contig_file)
        exit(1)
  
    checkv_cmd = [opt.checkv,
                  'end_to_end', opt.virsort2_pass1_contig_file,
                  opt.checkv_dir,
                  '-t', str(opt.num_cpu_threads),
                  '-d', opt.db_dir]

    run_sub_command(checkv_cmd, 'Run checkv')

    # combine result
    if not os.path.exists(opt.checkv_proviruses_contig):
        logger.error('can not find checkv result file: ' + opt.checkv_proviruses_contig)
        raise Usage('quit program')

    if not os.path.exists(opt.checkv_viruses_contig):
        logger.error('can not find checkv result file: ' + opt.checkv_viruses_contig)
        raise Usage('quit program')

    combine_cmd = ['cat', opt.checkv_proviruses_contig, opt.checkv_viruses_contig, ''.join(['>', opt.checkv_final_contig])]
    run_sub_command(combine_cmd, 'Combine checkv results')

@check_point
def viral_screen():
    cmd = [opt.viral_screen_script, 
           '-i', opt.input_contig_file,
           '-s', opt.virsort2_pass1_score_file,
           '-c', opt.checkv_contamination_table,
           '-o', opt.out_dir, '-p', opt.out_prefix]
    
    run_sub_command(cmd, "Viral screen based on virsorter and checkv results")


def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        start_time = time.time()
        setup_logger()
        parse_option(argv[1:])
        setup_output_dir()
    
        check_and_correct_option()

        check_input()
        opt.dump()

        logger.info('Start pipeline. Number of CPU threads %d ' % opt.num_cpu_threads)

        if opt.do_not_split:
            virsort_pass1()
        else:
            parallel_virsort_pass1()

        checkv()
        viral_screen()

        if not opt.keep_tmp_files and os.path.exists(opt.temp_dir):
            shutil.rmtree(opt.temp_dir)

        open(os.path.join(opt.out_dir, 'done'), 'w').close()

        logger.info('ALL DONE. Time elapsed: %f seconds ' % (time.time() - start_time))

    except Usage as usg:
        print(sys.argv[0].split('/')[-1] + ': ' + str(usg.msg), file=sys.stderr)
        exit(1)


if __name__ == '__main__':
    main()
