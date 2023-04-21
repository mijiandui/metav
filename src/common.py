
# -------------------------------------------------------------------------
# METAV
# Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics Limited
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

import os
import subprocess

def check_output(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    out, _ = p.communicate()
    out = out.rstrip().decode('utf-8')
    assert p.wait() == 0
    return out


def abspath(path):
    return os.path.abspath(os.path.expanduser(path))


def mkdir_if_not_exists(path):
    if not os.path.exists(path):
        os.mkdir(path)


def remove_if_exists(file_name):
    if os.path.exists(file_name):
        os.remove(file_name)


def detect_available_mem():
    try:
        psize = os.sysconf('SC_PAGE_SIZE')
        pcount = os.sysconf('SC_PHYS_PAGES')
        if psize < 0 or pcount < 0:
            raise SystemError
        return psize * pcount
    except ValueError:
        if sys.platform.find("darwin") != -1:
            return int(float(os.popen("sysctl hw.memsize").readlines()[0].split()[1]))
        elif sys.platform.find("linux") != -1:
            return int(float(os.popen("free").readlines()[1].split()[1]) * 1024)
        else:
            raise

def inpipe_cmd(file_name):
    if file_name.endswith('.gz'):
        return 'gzip -cd ' + file_name
    elif file_name.endswith('.bz2'):
        return 'bzip2 -cd ' + file_name
    else:
        return ''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


class EarlyTerminate(Exception):
    def __init__(self, kmer_size):
        self.kmer_size = kmer_size


class Checkpoint:

    def __init__(self, logger):
        self._current_checkpoint = 0
        self._logged_checkpoint = None
        self._file = None
        self.logger = logger

    def set_file(self, file):
        self._file = file

    def __call__(self, func):
        def checked_or_call(*args, **kwargs):
            if self._logged_checkpoint is None or self._current_checkpoint > self._logged_checkpoint:
                func(*args, **kwargs)
                with open(self._file, 'a') as cpf:
                    print(str(self._current_checkpoint) + '\t' + 'done', file=cpf)
            else:
                self.logger.info('passing check point {0}'.format(self._current_checkpoint))
            self._current_checkpoint += 1

        return checked_or_call

    def load_for_continue(self):
        self._logged_checkpoint = -1
        if os.path.exists(self._file):
            with open(self._file, 'r') as f:
                for line in f:
                    a = line.strip().split()
                    if len(a) == 2 and a[1] == 'done' and int(a[0]) > self._logged_checkpoint:
                        self._logged_checkpoint = int(a[0])


