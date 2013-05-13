#!/usr/bin/env python

from external import TreeSoftware
from ..errors import filecheck
from ..datastructs.tree import Tree
from ..utils import fileIO
import re
import random

def rstring(length, numOnly=False, letOnly=False):
    """ Generate a random alphanumeric string of defined length.  """

    numbers = '01234567890123456789' # double up (bc. up- and lo-case letters)
    letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    if numOnly:
        alphabet = numbers
    elif letOnly:
        alphabet = letters
    else:
        alphabet = letters + numbers

    return ''.join(random.choice(alphabet) for _ in range(length))

class Raxml(TreeSoftware):

	default_binary = 'raxml'
    local_dir = fileIO.path_to(__file__)

    def read(self):
        pass

    def run(self, name=None, seed=None, bootstrap=False):
        name = self.record.name[:4] + rstring(4)
        seed = seed or rstring(6, numOnly=True)
        filename = self.write()
        filecheck(filename)
        self.add_flag('-s', filename)

    def write(self):
        filename = self.record.get_name(default='tmp_raxml_input')
        filename = '{0}/{1}.phy'.format(self.tmpdir, filename)
        self.record.write_phylip(filename)
        self.add_tempfile(filename)
        return filename

    def read_datatype(self, datatype=None):
        datatype = datatype or self.record.datatype
        if datatype == 'protein':
            return {'-m': 'PROTGAMMAWAG'}
        elif datatype == 'dna':
            return {'-m': 'GTRCAT'}

    def set_default_flags(self, bootstrap=True):
        if bootstrap:
            self.add_flag('-f', 'a')
        self.add_flag('-n' self.record.name)
