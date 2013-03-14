#!/usr/bin/env python

from external import ExternalSoftware
from errors import filecheck, FileError
from utils import fileIO

class Darwin(ExternalSoftware):

    """ Run commands through Darwin language"""

    default_binary = 'darwin'

    def __init__(self, tmpdir='/tmp'):
        self.tmpdir = tmpdir
        self.comfile = '{0}/darcom.drw'.format(self.tmpdir)
        self.outfile = '{0}/output.drw'.format(self.tmpdir)

    def write(self, command):
        command = command.rstrip()
        writer = fileIO.fwriter(self.comfile)
        writer.write(command)
        if not command.endswith('\nquit;'):
            writer.write('\nquit;')
        writer.write('\n')
        writer.close()

    def execute(self, verbosity):
        comfile = filecheck('{0}/darcom.drw'.format(self.tmpdir))
        cmd = ('echo '
               '"outf := \'{0}\'; ReadProgram(\'{1}\');" '
               '| darwin'.format(self.outfile, comfile))
        if verbosity > 0:
            print cmd
        return fileIO.subprocess(cmd)

    def read(self):
        filecheck(self.outfile)
        reader = fileIO.freader(self.outfile)
        result = reader.read()
        reader.close()
        return result

    def clean(self):
        for f in (self.comfile, self.outfile):
            try:
                fileIO.delete(f)
            except FileError:
                continue

    def run(self, command, verbosity=0):
        self.write(command)
        (stdout, stderr) = self.execute(verbosity)
        if verbosity > 0:
            print stdout, stderr
        result = self.read()
        self.clean()
        return result
