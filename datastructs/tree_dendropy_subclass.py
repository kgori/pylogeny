#!/usr/bin/env python

import dendropy 
import re

class Tree(dendropy.Tree):
    """Augmented version of dendropy Tree class"""

    name_regex = re.compile('([A-Za-z0-9\-_]+).([A-Za-z0-9\-_]+)(?=_phyml_)')
    score_regex = re.compile('(?<=Log-likelihood: ).+')

    def __init__(
        self, 
        newick=None,
        score=0,
        output=None,
        program=None,
        name=None,
        ):
        super(Tree, self).__init__()
        if newick:
            self.read_from_string(newick, 'newick')
        self.score = score
        self.output = output
        self.program = program
        self.name = name    

    def __repr__(self):
        return '{0}{1}'.format(self.__class__.__name__,
                               (self.newick if self.newick else '(None)'))

    def __str__(self):
        """ Represents the object's information inside a newick comment, so is
        still interpretable by a (good) newick parser """

        s = '[Tree Object:\n'
        if self.name:
            s += 'Name:\t{0}\n'.format(self.name)

        s += ('Program:\t{0}\n'
              'Score:\t{1}\n'
              'Rooted:\t{2}\n'
              'Tree:\t]{3}\n'.format(self.program, self.score, 
                    self.rooted, self.newick))

        return s

    def __len__(self):
        return len(self.leaf_nodes())

    def __and__(self, other):
        return self.intersection(other)

    def copy(self):
        copy = self.__class__(self.newick, self.score, self.output, 
            self.program, self.name)
        return copy

    @property
    def newick(self):
        n = self.as_newick_string()
        if n:
            return (n if n.endswith(';') else n + ';')
        return n

    @property
    def rooted(self):
        return len(self.seed_node.child_nodes()) == 2 if self.newick else None

    @property
    def labels(self):
        return set([n.taxon.label for n in self.leaf_nodes()]) 

    @classmethod
    def bifurcate_base(cls, newick):
        t = cls(newick)
        t.resolve_polytomies()
        return t.newick

    @classmethod
    def trifurcate_base(cls, newick):
        t = cls(newick)
        t.deroot()
        return t.newick   

    def intersection(self, other):
        taxa1 = self.labels()
        taxa2 = other.labels()
        return taxa1 & taxa2

    def prune_to_subset(self, subset, inplace=False):
        if not subset.issubset(self.labels):
            print '"subset" is not a subset'
            return
        if not inplace:
            t = self.copy()
        else:
            t = self
        t.retain_taxa_with_labels(subset)
        return t

    def pruned_pair(self, other, inplace=False):
        assert isinstance(other, self.__class__)


    def scale(self, factor, inplace=True):
        if not inplace:
            t = self.copy()
        else:
            t = self
        t.scale_edges(factor)
        return t

    def write_to_file(
        self,
        outfile,
        metadata=False,
        suppress_NHX=False,
        ):
        """ Writes a string representation of the object's contents to file.
        This can be read using read_from_file to reconstruct the Tree object, if
        metadata is included (i.e. metadata=True) """

        with open(outfile, 'w') as writer:
            if metadata:
                writer.write(str(self))
            else:
                writeable = self.newick
                if suppress_NHX:
                    if writeable.startswith('[&R] '):
                        writeable = writeable[5:]
                if not writeable.endswith('\n'):
                    writeable += '\n'
                writer.write(writeable)        
        return outfile

