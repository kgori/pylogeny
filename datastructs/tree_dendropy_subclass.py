#!/usr/bin/env python

import dendropy
import re
import random
from ..errors import FileError, filecheck


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
        **kwargs
        ):

        super(Tree, self).__init__()
        if newick:
            self.read_from_string(newick, 'newick', **kwargs)
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

        s += \
            '''Program:\t{0}
Score:\t{1}
Rooted:\t{2}
Tree:\t]{3}
'''.format(self.program,
                self.score, self.rooted, self.newick)

        return s

    def __len__(self):
        """ Number of leaves on the Tree. For total branch length use 
        self.length()"""
        return len(self.leaf_nodes())

    def __and__(self, other):
        """ Overloads & operator:

        'self & other' is equivalent to 'self.intersection(other)''
        """
        return self.intersection(other)

    def copy(self):
        """ Returns an independent copy of self """
        copy = self.__class__(self.newick, self.score, self.output,
                              self.program, self.name)
        return copy

    @property
    def newick(self):
        n = self.as_newick_string()
        if n:
            return (n if n.endswith(';') else n + ';')
        return n

    @newick.setter
    def newick(self, newick_string):
        if self.newick:
            print 'Newick string already loaded: {0}'.format(self.newick)
            return
        self.read_from_string(newick_string, 'newick')

    @property
    def rooted(self):
        """ Predicate testing for rootedness by checking for a bifurcation 
        at the root. """
        return (len(self.seed_node.child_nodes()) == 2 if self.newick else None)

    @property
    def labels(self):
        """ Returns the taxon set of the tree (same as the label- or 
        leaf-set) """
        return set([n.taxon.label for n in self.leaf_nodes()])

    @classmethod
    def bifurcate_base(cls, newick):
        """ Rewrites a newick string so that the base is a bifurcation
        (rooted tree) """
        t = cls(newick)
        t.resolve_polytomies()
        return t.newick

    @classmethod
    def trifurcate_base(cls, newick):
        """ Rewrites a newick string so that the base is a trifurcation
        (usually means an unrooted tree) """
        t = cls(newick)
        t.deroot()
        return t.newick

    def intersection(self, other):
        """ Returns the intersection of the taxon sets of two Trees """
        taxa1 = self.labels()
        taxa2 = other.labels()
        return taxa1 & taxa2

    def prune_to_subset(self, subset, inplace=False):
        """ Prunes the Tree to just the taxon set given in `subset` """
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
        """ Returns two trees pruned to the intersection of their taxon sets """
        assert isinstance(other, self.__class__)

    def scale(self, factor, inplace=True):
        """ Multiplies all branch lengths by factor. """
        if not inplace:
            t = self.copy()
        else:
            t = self
        t.scale_edges(factor)
        return t

    def strip(self, inplace=False):
        """ Sets all edge lengths to None """
        if not inplace:
            t = self.copy()
        else:
            t = self
        for e in t.preorder_edge_iter():
            e.length = None
        return t        

    def randomise_branch_lengths(
        self,
        i=(1, 1),
        l=(1, 1),
        distribution_func=random.gammavariate,
        inplace=False,
        ):
        """ Replaces branch lengths with values drawn from the specified
        distribution_func. Parameters of the distribution are given in the 
        tuples i and l, for interior and leaf nodes respectively. """

        if not inplace:
            t = self.copy()
        else:
            t = self

        for n in t.preorder_node_iter():
            if n.is_internal():
                n.edge.length = max(0, distribution_func(*i))
            else:
                n.edge.length = max(0, distribution_func(*l))
        return t

    def get_inner_edges(self):
        """ Returns a list of the internal edges of the tree. """
        inner_edges = [e for e in self.preorder_edge_iter() if e.is_internal()
                       and e.head_node and e.tail_node]
        return inner_edges

    def get_children(tree, inner_edge):
        """ Given an edge in the tree, returns the child nodes of the head and
        the tail nodes of the edge, for instance:

            A      C    | A, B, C and D are the children of the edge --->,
             \    /     | C and D are the head node children, and A and B
              --->      | are the tail node children.   
             /    \     
            B      D    | Output: {'head': [<C>, <D>], 'tail': [<A>, <B>]}

        N1: Edges are directional in dendropy trees. The head node of an
        edge is automatically a child of the tail node, but we don't want this.
        """

        h = inner_edge.head_node
        t = inner_edge.tail_node
        if not tree.seed_node == t:
            original_seed = tree.seed_node
            tree.reseed_at(t)
        else:
            original_seed = None
        head_children = h.child_nodes()
        tail_children = list(set(t.child_nodes()) - set([h]))  # See N1
        if original_seed:
            tree.reseed_at(original_seed)
        return {'head': head_children, 'tail': tail_children}

    def nni(
        self,
        edge,
        head_subtree,
        tail_subtree,
        ):
        """ Nearest-neighbour interchange (NNI) operation.
        
        An edge in the tree has two or more subtrees at each end (ends are
        designated 'head' and 'tail'). The NNI operation exchanges one of the
        head subtrees for one of the tail subtrees, as follows:

            A      C                        C      A    | Subtree A is exchanged
             \    /        +NNI(A,C)         \    /     | with subtree C.
              --->        ==========>         --->      |
             /    \                          /    \     |
            B      D                        B      D
        
        This operation acts on the tree inplace.
        """

        original_seed = self.seed_node
        head = edge.head_node
        tail = edge.tail_node
        self.reseed_at(tail)
        assert head_subtree.parent_node == head
        assert tail_subtree.parent_node == tail
        head.remove_child(head_subtree)
        tail.remove_child(tail_subtree)
        head.add_child(tail_subtree)
        tail.add_child(head_subtree)
        self.reseed_at(original_seed)
        self.update_splits()

    def rnni(self, inplace=False):
        """ Applies a NNI operation on a randomly chosen edge. If called with
        inplace=True this will alter the structure of the calling Tree. """
        if inplace:
            tree = self
        else:
            tree = self.copy()

        e = random.choice(tree.get_inner_edges())
        children = tree.get_children(e)
        (h, t) = (random.choice(children['head']), random.choice(children['tail'
                  ]))
        tree.nni(e, h, t)
        return tree

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

    def ntaxa(self):
        return len(self)

    @classmethod
    def new_tree_from_phyml_results(
        cls,
        tree_file,
        stats_file,
        name=None,
        program='phyml',
        ):
        """ Given the usual phyml output files - xxx_phyml_tree.txt and 
        xxx_phyml_stats.txt, instantiates a Tree from the information
        in the phyml files. """
        #newick score output program name
        
        exit = False
        for f in (tree_file, stats_file):
            try:
                filecheck(f)
            except FileError, e:
                print e
                exit = True

        if exit:
            print 'Results were not loaded'
            raise FileError()
        
        if not name:
            name = cls.name_regex.search(tree_file).group(1)
        newick = open(tree_file).read()
        stats = open(stats_file).read()
        score = cls.score_regex.search(stats).group(0)
        score = float(score) if score else None

        return cls(newick=newick, score=score, output=stats, program=program,
            name=name)

    def _unify_taxon_sets(self, other):
        if other.taxon_set is not self.taxon_set:
            return self.__class__(other.newick, taxon_set=self.taxon_set)
        else:
            return other

    def rfdist(self, other):
        cp = self._unify_taxon_sets(other)
        return self.symmetric_difference(cp)

    def eucdist(self, other):
        cp = self._unify_taxon_sets(other)
        return self.euclidean_distance(cp)

    def wrfdist(self, other):
        cp = self._unify_taxon_sets(other)
        return self.robinson_foulds_distance(cp)
