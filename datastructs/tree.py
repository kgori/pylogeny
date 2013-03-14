#!/usr/bin/env python

from ..utils import fileIO


class Tree(object):

    def __init__(
        self,
        newick,
        score=0,
        output='',
        program='',
        name='tree_object',
        ):

        self.newick = newick
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
              'Tree:\t]{2}\n'.format(self.program, self.score, self.newick))

        return s

    def __eq__(self, other):

        return all([self.name == other.name, self.newick == other.newick,
                   self.program == other.program, self.score == other.score,
                   self.output == other.output])

    def copy(self):
        copy = self.__new__(type(self))
        copy.__dict__ = {key: value for (key, value) in self.__dict__.items()}
        return copy

    def write_to_file(
        self,
        outfile,
        metadata=False,
        suppress_NHX=False,
        ):
        """ Writes a string representation of the object's contents to file.
        This can be read using read_from_file to reconstruct the Tree object, if
        metadata is included (i.e. metadata=True) """

        with fileIO.fwriter(outfile) as writer:

            if metadata:
                writer.write(str(self))
            else:

                writeable = self.newick
                if suppress_NHX:
                    writeable.lstrip('[&R] ')
                if not writeable.endswith('\n'):
                    writeable += '\n'
                writer.write(writeable)

        return outfile
