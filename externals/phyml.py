#!/usr/bin/env python

from external import ExternalSoftware, TreeSoftware
from ..errors import filecheck, optioncheck
from ..datastructs.tree import Tree
from ..utils import fileIO
from ..utils.printing import print_and_return
from bsub import bsub
import os
import re
import shutil
import tempfile


ANALYSES = ['tlr', 'lr', 'l', 'r', 'ml', 'full', 'nj', 'bionj', 'bionj+', 'lk']

class LSFPhyml(ExternalSoftware):

    default_binary = 'phyml'
    score_regex = re.compile('(?<=Log-likelihood: ).+')
    local_dir = fileIO.path_to(__file__)

    def __init__(self, records, tmpdir, supplied_binary=''):
        super(LSFPhyml, self).__init__(tmpdir, supplied_binary)
        self.records = records
        self.temp_dirs = self.setup_temp_dirs()
        self.phyml_objects = self.setup_phyml_objects()

    @property
    def records(self):
        return self._records

    @records.setter
    def records(self, records):
        self._records = records

    def write(self):
        pass

    def setup_temp_dirs(self):
        temp_dirs = [tempfile.mkdtemp(dir=self.tmpdir) for rec in self.records]
        return temp_dirs

    def setup_phyml_objects(self):
        phyml_objects = [Phyml(rec, td, self.binary)
                         for (rec, td) in zip(self.records, self.temp_dirs)]
        return phyml_objects

    def get_command_strings(self, analysis='ml'):
        return [phyml.run(analysis, dry_run=True)
                for phyml in self.phyml_objects]

    def launch_lsf(self, command_strings, verbose=False, output='/dev/null'):
        curr_dir = os.getcwd()
        os.chdir(self.tmpdir)
        job_ids = [bsub('phyml_task',
                        o='/dev/null',
                        e='/dev/null',
                        verbose=verbose)(cmd).job_id
                   for cmd in command_strings]
        bsub.poll(job_ids)
        os.chdir(curr_dir)

    def read(self, analysis):
        self.trees = []
        for phyml in self.phyml_objects:
            (tree, stats) = phyml.read()
            try:
                score = float(self.score_regex.search(stats).group(0))
            except:
                score = 0
            tree_object = Tree(newick=tree, score=score, program=analysis,
                               name=phyml.record.name, output=stats)
            self.trees.append(tree_object)
        return self.trees

    def clean(self):
        for phyml in self.phyml_objects:
            phyml.clean()
        for d in self.temp_dirs:
            shutil.rmtree(d)
        for f in os.listdir(self.tmpdir):
            if (f.endswith('.out') or f.endswith('.err')):
                os.remove(f)

    def run(self, analysis, verbose=False):
        command_strings = self.get_command_strings(analysis)
        self.launch_lsf(command_strings, verbose)
        trees = self.read(analysis)
        if len(trees) == len(self.records):
            self.clean()
        return trees

    def call(self):
        pass


class Phyml(TreeSoftware):

    """ __init__ takes a Seq sequence record as
    first positional argument, tmpdir as second, and supplied_binary=
    as keyword argument """

    default_binary = 'phyml'
    score_regex = re.compile('(?<=Log-likelihood: ).+')
    local_dir = fileIO.path_to(__file__)

    def read(self, filename=None):
        filename = filename or self.filename
        tree_filename = filename + '_phyml_tree.txt'
        stats_filename = filename + '_phyml_stats.txt'
        self.add_tempfile(tree_filename)
        self.add_tempfile(stats_filename)
        with open(tree_filename) as treefile:
            with open(stats_filename) as statsfile:
                return (treefile.read(), statsfile.read())

    def run(self, analysis=None, verbosity=0,
            **kwargs):
        if analysis:
            self.set_default_flags(analysis)
        else:
            analysis = fileIO.basename(self.binary)
        optioncheck(analysis, ANALYSES)
        if verbosity > 1:
            print self.flags
            print 'Writing tempfiles to', self.tmpdir
        filename = self.write()
        filecheck(filename)
        self.add_flag('-i', filename)

        # OUTPUT TO USER
        if verbosity == 1:
            print_and_return('Running phyml on {0}'.format(self.record.name))
        elif verbosity > 1:
            print 'Running phyml on {0}'.format(self.record.name)

        # DRY RUN - just get command string
        if kwargs.get('dry_run', False):
            cmd = self.call(verbose=(True if verbosity > 1 else False),
                dry_run=True)
            return cmd

        # RUN PHYML
        (stdout, stderr) = self.call(verbose=(True if verbosity > 1 else False))
        (tree, stats) = self.read(filename)
        try:
            score = float(self.score_regex.search(stats).group(0))
        except:
            print tree
            print stats
        if verbosity > 1:
            print 'Cleaning tempfiles'
        self.clean()
        tree_object = Tree(newick=tree, score=score, program=analysis,
            name=self.record.name, output=stats, **kwargs)
        if kwargs.get('set_as_record_tree', True):
            self.record.tree = tree_object
        if verbosity > 1:
            print 'Done.'
        return tree_object

    def write(self):
        record_name = self.record.get_name(default='tmp_phyml_input')
        with tempfile.NamedTemporaryFile(prefix=self.record.name,
                                         suffix='.phy',
                                         dir=self.tmpdir,
                                         delete=False) as file_:
            filename = file_.name
            file_.write(self.record.write_phylip('pipe'))
        self.filename = filename
        # filename = '{0}/{1}.phy'.format(self.tmpdir, record_name)
        self.add_tempfile(filename)
        return filename

    def read_datatype(self, datatype=None):
        datatype = datatype or self.record.datatype
        if datatype == 'protein':
            return {'-d': 'aa', '-m': 'WAG'}
        elif datatype == 'dna':
            return {'-d': 'nt', '-m': 'GTR'}

    def set_default_flags(self, analysis='ml', datatype=None):

        defaults = self.read_datatype(datatype=datatype)
        if defaults:
            defaults['-a'] = 'e'
            defaults['-b'] = 0
            defaults['-c'] = 4
            defaults['-q'] = ''
            defaults['--no_memory_check'] = ''
            defaults['--quiet'] = ''
            if analysis == 'ml' or analysis == 'full' or analysis == 'tlr':
                defaults['-o'] = 'tlr'
            elif analysis == 'nj' or analysis == 'bionj':
                defaults['-o'] = 'n'
            elif analysis == 'lr' or analysis == 'bionj+':
                defaults['-o'] = 'lr'
            elif analysis == 'l':
                defaults['-o'] = 'l'
            elif analysis == 'r':
                defaults['-o'] = 'r'
            elif analysis == 'lk':
                defaults['-o'] = 'n'

            for flag in defaults:
                self.add_flag(flag, defaults[flag])

def runPhyml(rec, tmpdir, analysis, verbosity=0, tree=None, **kwargs):
    optioncheck(analysis, ANALYSES)
    p = Phyml(rec, tmpdir)
    if analysis == 'lk' and tree is not None:
        tree_name = (tree.name if tree.name else 'tmp_tree')
        tmp_treefile = '{0}/{1}.nwk'.format(tmpdir, tree_name)
        tree.write_to_file(tmp_treefile)
        p.add_tempfile(filecheck(tmp_treefile))
        p.add_flag('-u', tmp_treefile)
    return p.run(analysis, verbosity, **kwargs)

def runLSFPhyml(records, tmpdir, analysis, verbosity, **kwargs):
    optioncheck(analysis, ANALYSES)
    lsfphyml = LSFPhyml(records, tmpdir)
    trees = lsfphyml.run(analysis, verbose=(True if verbosity > 0 else False))
    return trees
