#!/usr/bin/env python

import dendropy 

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
        copy = self.__new__(type(self))
        copy.__dict__ = {key: value for (key, value) in self.__dict__.items()}
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

    def labels(self):
        return set([n.taxon.label for n in self.leaf_nodes()])    

    def intersection(self, other):
        taxa1 = self.labels()
        taxa2 = other.labels()
        return taxa1 & taxa2

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

