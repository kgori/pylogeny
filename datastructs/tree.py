#!/usr/bin/env python

import dendropy
import re
import random
from ..errors import FileError, filecheck

class TreeEdgeError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg

def rootcheck(edge, msg='This is the root edge'):
    """ Raises error if edge is the root edge (has no tail node) """
    if not edge.tail_node:
        raise TreeEdgeError(msg)

def edge_length_check(length, edge):
    """ Raises error if length is not in interval [0, edge.length] """
    try:
        assert 0 <= length <= edge.length
    except AssertionError:
        if length < 0:
            raise TreeEdgeError('Negative edge-lengths are disallowed')
        raise TreeEdgeError(
            'This edge isn\'t long enough to prune at length {0}\n'
            '(Edge length = {1})'.format(length, edge.length))

def cast(dendropy_tree):
    """ Cast dendropy.Tree instance as Tree instance """
    return Tree(dendropy_tree.as_newick_string() + ';')

class LGT(object):

    """ This class provides a random LGT implementation using SPRs from an older
    version of this software. I've kept it for now because it works, and to 
    reimplement from scratch would be inconvenient, but is most definitely on 
    the 'to do' list """ 

    def __init__(self, tree, verbosity=0):
        self.tree = tree
        self.verbosity = verbosity

    def get_blocks(self):
        dists = []
        blocks = {}
        rootlen = self.tree.seed_node.edge_length
        
        for n in self.tree.preorder_node_iter():
            node_height = n.distance_from_root() - rootlen
            if not n.parent_node:
                continue
                root_height = n.distance_from_root()
                self.tree_height = root_height + n.distance_from_tip()
                parent_height = 0
            else:
                parent_height = n.parent_node.distance_from_root() - rootlen
            node_height = round(node_height, 8)
            parent_height = round(parent_height, 8)

            if not node_height in blocks:
                blocks[node_height] = []

            dists.append((n, parent_height, node_height))

        for time in blocks:
            for (node, parent_h, node_h) in dists:
                if parent_h < time <= node_h:
                    blocks[time].append(node)

        dists.sort(key=lambda x: x[2])
        if self.verbosity > 0:
            print blocks, dists
        return (blocks, dists)

    def weight_by_branches(self, blocks):
        intervals = sorted(blocks.keys())
        weighted_intervals = [0] + [None] * (len(intervals) - 1)
        for i in range(1, len(intervals)):
            time_range = intervals[i] - intervals[i - 1]
            num_branches = len(blocks[intervals[i]])
            weighted_range = time_range * num_branches
            weighted_intervals[i] = weighted_range + weighted_intervals[i
                    - 1]
        if self.verbosity > 0:
            print weighted_intervals
        return weighted_intervals

    def get_time(self, blocks, weights=None):
        d = sorted(blocks.keys())
        if weights:
            samp = random.uniform(weights[0], weights[-1])
            for i in range(len(weights) - 1):
                if weights[i + 1] >= samp > weights[i]:
                    interval = weights[i + 1] - weights[i]
                    proportion = (samp - weights[i]) / interval
                    break
            drange = d[i + 1] - d[i]
            time = drange * proportion + d[i]
        else:
            time = random.uniform(d[0], d[-1])

        if self.verbosity > 0: 
            print 'LGT event at time: {0}'.format(time)

        return time

    def choose_prune_and_regraft_nodes(self, time, dists):
        matching_branches = [x for x in dists if x[1] < time < x[2]]
        if self.verbosity > 0:
            print matching_branches

        prune = random.sample(matching_branches, 1)[0]
        if self.verbosity > 0:
            print prune
        siblings = prune[0].sister_nodes()
        for br in matching_branches:
            if br[0] in siblings:
                matching_branches.remove(br)

        matching_branches.remove(prune)
        if self.verbosity > 0:
            print matching_branches

        if matching_branches == []:
            if self.verbosity > 0:
                print 'No non-sibling branches available'
            return (None, None)

        regraft = random.sample(matching_branches, 1)[0]

        prune_taxa = [n.taxon.label for n in prune[0].leaf_iter()]
        regraft_taxa = [n.taxon.label for n in regraft[0].leaf_iter()]
        if self.verbosity:
            print 'Donor group = {0}'.format(regraft_taxa)
            print 'Receiver group = {0}'.format(prune_taxa)
        return (prune, regraft)  

    def lgt(self, time=None):
        """ Known bugs: will hang at P1 if time in range [0 - depth of 1st node]
        """

        (blocks, dists) = self.get_blocks()
        if not time:
            weights = self.weight_by_branches(blocks)
            time = self.get_time(blocks, weights)
        (p, r) = (None, None)
        while (p, r) == (None, None): #P1
            (p, r) = self.choose_prune_and_regraft_nodes(time, dists)
        pruning_edge = p[0].edge
        regrafting_edge = r[0].edge
        length1 = p[2] - time 
        length2 = r[2] - time
        return pruning_edge, regrafting_edge, length1, length2

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
            if self.rooted:
                self.is_rooted = True
                self.update_splits()
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

        s += ('Program:\t{0}'
            'Score:\t{1}'
            'Rooted:\t{2}'
            'Tree:\t]{3})'.format(self.program,
                self.score, self.rooted, self.newick))

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

    def randomise_labels(
        self,
        inplace=False,
        ):
        """ Shuffles the leaf labels, but doesn't alter the tree structure """

        if not inplace:
            t = self.copy()
        else:
            t = self

        names = t.labels
        random.shuffle(list(names))
        for l in t.leaf_iter():
            l.taxon_label = names.pop()
        return t

    def get_inner_edges(self):
        """ Returns a list of the internal edges of the tree. """
        inner_edges = [e for e in self.preorder_edge_iter() if e.is_internal()
                       and e.head_node and e.tail_node]
        return inner_edges

    def get_nonroot_edges(self):
        return [e for e in self.preorder_edge_iter()
             if e.head_node and e.tail_node]

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
        """ *Inplace* Nearest-neighbour interchange (NNI) operation.
        
        An edge in the tree has two or more subtrees at each end (ends are
        designated 'head' and 'tail'). The NNI operation exchanges one of the
        head subtrees for one of the tail subtrees, as follows:

            A      C                        C      A    | Subtree A is exchanged
             \    /        +NNI(A,C)         \    /     | with subtree C.
              --->        ==========>         --->      |
             /    \                          /    \     |
            B      D                        B      D
        
        
        """
        
        # This implementation works on unrooted Trees. If the input Tree is
        # rooted, it is derooted in a way that allows the root to be
        # re-established, using the reversible_deroot() method
        reroot = False                 
        if self.rooted:                 
            reroot = True
            rooting_data = self.reversible_deroot() 

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
        if reroot:
            self.reroot_at_edge(*rooting_data)
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

    def prune(self, edge, length=None):
        """ Prunes a subtree from the main Tree, retaining an edge length
        specified by length (defaults to entire length). The length is sanity-
        checked by edge_length_check, to ensure it is within the bounds 
        [0, edge.length].

        Returns the basal node of the pruned subtree. """

        length = length or edge.length
        edge_length_check(length, edge)

        n = edge.head_node
        self.prune_subtree(n)
        n.edge_length = length
        return n

    def regraft(self, edge, node, length=None):
        """ Grafts a node onto an edge of the Tree, at a point specified by
        length (defaults to middle of edge). """

        length = length or edge.length/2. # Length measured from head to tail
        edge_length_check(length, edge)
        rootcheck(edge, 'SPR regraft is not allowed on the root edge')

        t = edge.tail_node
        h = edge.head_node
        new = t.new_child(edge_length=edge.length-length)
        t.remove_child(h)
        new.add_child(h, edge_length=length)
        new.add_child(node)
        self.update_splits()

    def spr(self, pruning_edge, regrafting_edge, length1=None, length2=None):
        """ Combines self.prune() and self.regraft() methods to perform
        Subtree-Prune and Regraft (SPR) operation. 

        Notes
        N1: Sanity check that we don't try to regraft to pruned subtree onto
        an edge that is in the pruned subtree
        N2: Adjustment to account for the fact the the pruned edge's parent
        edge is deleted in the prune operation.
        """
        try: # See note N1
            descendents = [n.edge for n in 
                pruning_edge.head_node.preorder_iter()]
            assert regrafting_edge not in descendents
        except AssertionError:
            raise TreeEdgeError('Regraft edge is on the pruned subtree')

        sister_nodes = pruning_edge.head_node.sister_nodes() # See note N2
        if (regrafting_edge == pruning_edge.tail_node.edge
            and len(sister_nodes) == 1):
            regrafting_edge = sister_nodes[0].edge
            length2 += sister_nodes[0].edge_length # end note N2
        
        n = self.prune(pruning_edge, length1)        
        self.regraft(regrafting_edge, n, length2)

    def rspr(self, inplace=False, disallow_sibling_sprs=False):
        """ Random SPR, with prune and regraft edges chosen randomly, and
        lengths drawn uniformly from the available edge lengths.

        N1: disallow_sibling_sprs prevents sprs that don't alter the topology
        of the tree """
        if not inplace:
            tree = self.copy()
        else:
            tree = self
        
        edges = [e for e in tree.preorder_edge_iter()
                 if e.head_node and e.tail_node]
        pruning_edge = random.choice(edges)
        
        if disallow_sibling_sprs: # See note N1
            for e in pruning_edge.adjacent_edges:
                edges.remove(e)

        for n in pruning_edge.head_node.preorder_iter():
            edges.remove(n.edge)

        regrafting_edge = random.choice(edges)

        length1 = random.uniform(0, pruning_edge.length)
        length2 = random.uniform(0, regrafting_edge.length)

        tree.spr(pruning_edge, regrafting_edge, length1, length2)
        return tree

    def reversible_deroot(self):
        """ Stores info required to restore rootedness to derooted Tree. Returns
        the edge that was originally rooted, the length of e1, and the length 
        of e2.

        Dendropy Derooting Process:
        In a rooted tree the root node is bifurcating. Derooting makes it 
        trifurcating. 

        Call the two edges leading out of the root node e1 and e2.
        Derooting with Tree.deroot() deletes one of e1 and e2 (let's say e2), 
        and stretches the other to the sum of their lengths. Call this e3.

        Rooted tree:                   Derooted tree:
                 A                         A   B
                 |_ B                       \ /   
                /                            |
               /e1                           |e3 (length = e1+e2; e2 is deleted)
        Root--o               ===>           |
               \e2                     Root--o _ C   
                \ _ C                        |
                 |                           D
                 D

        Reverse this with Tree.reroot_at_edge(edge, length1, length2, ...)
        """
        root_edge = self.seed_node.edge
        lengths = dict([(edge, edge.length) for edge 
            in self.seed_node.incident_edges() if edge is not root_edge])
        self.deroot()
        reroot_edge = (set(self.seed_node.incident_edges()) 
                            & set(lengths.keys())).pop()
        return (reroot_edge, reroot_edge.length - lengths[reroot_edge], 
            lengths[reroot_edge])

    def rlgt(self, inplace=False, time=None, verbosity=0):
        """ Uses class LGT to perform random lateral gene transfer on 
        ultrametric tree """

        if not inplace:
            t = self.copy()
        else:
            t = self

        try:
            self.calc_node_ages() # Fails with ValueError if tree is not 
        except ValueError:        # ultrametric
            raise

        lgt = LGT(t)
        pr_ed, rg_ed, l1, l2 = lgt.lgt(time)
        t.spr(pr_ed, rg_ed, l1, l2)
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

    def ntaxa(self):
        return len(self)

    def _name_things(self):
        """ Easy names for debugging """
        edges = {}
        nodes = {None: 'root'}
        for n in self.postorder_node_iter():
            nodes[n] = '.'.join([str(x.taxon) for x in n.leaf_nodes()])
        for e in self.preorder_edge_iter():
            edges[e] = ' ---> '.join([nodes[e.tail_node], nodes[e.head_node]])

        r_edges = {value: key for key,value in edges.items()}
        r_nodes = {value: key for key,value in nodes.items()}
        return edges, nodes, r_edges, r_nodes

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

    ### DISTANCE CALCULATIONS
    # _unify_taxon_sets() might not be necessary - can't seem to reproduce
    # the bug that made me introduce it

    def _unify_taxon_sets(self, other):
        if other.taxon_set is not self.taxon_set:
            return self.__class__(other.newick, taxon_set=self.taxon_set)
        else:
            return other

    def rfdist_(self, other):
        cp = self._unify_taxon_sets(other)
        return self.symmetric_difference(cp)

    def rfdist(self, other):
        try:
            return self.symmetric_difference(other)
        except:
            print 'ERROR!'
            return self.rfdist_(other)

    def eucdist_(self, other):
        cp = self._unify_taxon_sets(other)
        return self.euclidean_distance(cp)

    def eucdist(self, other):
        try:
            return self.euclidean_distance(other)
        except:
            print 'ERROR!'
            return self.eucdist_(other)

    def wrfdist_(self, other):
        cp = self._unify_taxon_sets(other)
        return self.robinson_foulds_distance(cp)

    def wrfdist(self, other):
        try:
            return self.robinson_foulds_distance(other)
        except:
            print 'ERROR!'
            return self.wrfdist_(other)

    @classmethod
    def new_rtree(cls, nspecies=16, **kwargs):
        tg = TreeGen(nspecies, **kwargs)
        return tg.rtree()

    @classmethod
    def new_coal(cls, nspecies=16, **kwargs):
        tg = TreeGen(nspecies, **kwargs)
        return tg.coal()

    @classmethod
    def new_yule(cls, nspecies=16, **kwargs):
        tg = TreeGen(nspecies, **kwargs)
        return tg.yule()

    def sample_gene_tree(self, **kwargs):
        tg = TreeGen(template=self)
        return tg.gene_tree(**kwargs)

class TreeGen(object):

    def __init__(
        self,
        nspecies=16,
        names=None,
        template=None,
        cf=False,
        ):

        """ Generates a new Tree using a coalescent process (coal method), a
        Yule pure-birth process (yule method), a random tree (rtree), or by
        sampling a gene tree from a template species tree using a constrained
        Kingman coalescent.

        nspecies = number of taxa in the tree
        names = a list of leaf names (names will be generated if not supplied)
        template = a template species tree for drawing gene trees
        cf = set to true to generate leaf names from the list of character names
        from Cannon Fodder """

        self.nspecies = nspecies
        if cf:
            self.names = random.sample(cfnames, nspecies)
        else:
            self.names = names or ['Sp{0}'.format(i) for i in range(1, nspecies
                                   + 1)]
        if template and not isinstance(template, Tree):
            raise TypeError('template should be \'Tree\' object. Got',
                            type(tree))
        self.template = template
    
    def coal(self):
        taxon_set = dendropy.TaxonSet(self.names)
        return cast(dendropy.treesim.pure_kingman(taxon_set))

    def gene_tree(
        self,
        scale_to=None,
        population_size=1,
        trim_names=True,
        ):
        """ Using the current tree object as a species tree, generate a gene
        tree using the constrained Kingman coalescent process from dendropy. The
        species tree should probably be a valid, ultrametric tree, generated by
        some pure birth, birth-death or coalescent process, but no checks are
        made. Optional kwargs are: -- scale_to, which is a floating point value
        to scale the total tree tip-to-root length to, -- population_size, which
        is a floating point value which all branch lengths will be divided by to
        convert them to coalescent units, and -- trim_names, boolean, defaults
        to true, trims off the number which dendropy appends to the sequence
        name """

        t = self.template or self.yule()

        for leaf in tree.leaf_iter():
            leaf.num_genes = 1

        tree_height = tree.seed_node.distance_from_root() \
            + tree.seed_node.distance_from_tip()

        if scale_to:
            population_size = tree_height / scale_to

        for edge in tree.preorder_edge_iter():
            edge.pop_size = population_size

        gene_tree = dendropy.treesim.constrained_kingman(tree)[0]

        if trim_names:
            for leaf in gene_tree.leaf_iter():
                leaf.taxon.label = leaf.taxon.label.replace('\'', '').split('_'
                        )[0]

        return {'gene_tree': cast(gene_tree), 'species_tree': t}

    def rtree(self):
        m = self.yule()
        m.randomise_labels()
        return m.randomise_branch_lengths()

    def yule(self):
        taxon_set = dendropy.TaxonSet(self.names)
        return cast(dendropy.treesim.uniform_pure_birth(taxon_set))


cfnames = [
    'Jools', 'Jops', 'Stoo', 'Rj', 'Ubik', 'Cj', 'Chris', 'Pete',
    'Tadger', 'Hector', 'Elroy', 'Softy', 'Mac', 'Bomber', 'Stan', 'Tosh',
    'Brains', 'Norm', 'Buster', 'Spike', 'Browny', 'Murphy', 'Killer', 'Abdul',
    'Spotty', 'Goofy', 'Donald', 'Windy', 'Nifta', 'Denzil', 'Cedric', 'Alf',
    'Marty', 'Cecil', 'Wally', 'Pervy', 'Jason', 'Roy', 'Peewee', 'Arnie',
    'Lofty', 'Tubby', 'Porky', 'Norris', 'Bugsy', 'Greg', 'Gus', 'Ginger',
    'Eddy', 'Steve', 'Hugo', 'Zippy', 'Sonny', 'Willy', 'Mario', 'Luigi',
    'Bo', 'Johan', 'Colin', 'Queeny', 'Morgan', 'Reg', 'Peter', 'Brett',
    'Matt', 'Vic', 'Hut', 'Bud', 'Brad', 'Ashley', 'Les', 'Rex',
    'Louis', 'Pedro', 'Marco', 'Leon', 'Ali', 'Tyson', 'Tiger', 'Frank',
    'Reuben', 'Leyton', 'Josh', 'Judas', 'Aj', 'Lex', 'Butch', 'Bison',
    'Gary', 'Luther', 'Kermit', 'Brian', 'Ray', 'Freak', 'Leroy', 'Lee',
    'Banjo', 'Beaker', 'Basil', 'Bonzo', 'Kelvin', 'Ronnie', 'Rupert', 'Roo',
    'Dan', 'Jimmy', 'Bob', 'Don', 'Tommy', 'Eddie', 'Ozzy', 'Paddy',
    'Arnold', 'Tony', 'Teddy', 'Dom', 'Theo', 'Martin', 'Chunky', 'Jon',
    'Ben', 'Girly', 'Julian', 'Pizza', 'Ciaran', 'Jock', 'Gravy', 'Trendy',
    'Neil', 'Derek', 'Ed', 'Biff', 'Paul', 'Stuart', 'Randy', 'Loreta',
    'Suzie', 'Pumpy', 'Urmer', 'Roger', 'Pussy', 'Meat', 'Beefy', 'Harry',
    'Tiny', 'Howard', 'Morris', 'Thor', 'Rev', 'Duke', 'Micky', 'Chas',
    'Melony', 'Craig', 'Sidney', 'Parson', 'Rowan', 'Smelly', 'Dok', 'Stew',
    'Adrian', 'Pat', 'Iceman', 'Goose', 'Dippy', 'Viv', 'Fags', 'Bunty',
    'Noel', 'Bono', 'Edge', 'Robbie', 'Sean', 'Miles', 'Jimi', 'Gordon',
    'Val', 'Hobo', 'Fungus', 'Toilet', 'Lampy', 'Marcus', 'Pele', 'Hubert',
    'James', 'Tim', 'Saul', 'Andy', 'Silky', 'Simon', 'Handy', 'Sid',
    'George', 'Joff', 'Barry', 'Dick', 'Gil', 'Nick', 'Ted', 'Phil',
    'Woody', 'Wynn', 'Alan', 'Pip', 'Mickey', 'Justin', 'Karl', 'Maddog',
    'Horace', 'Harold', 'Gazza', 'Spiv', 'Foxy', 'Ned', 'Bazil', 'Oliver',
    'Rett', 'Scot', 'Darren', 'Moses', 'Noah', 'Seth', 'Buddah', 'Mary',
    'Pilot', 'Mcbeth', 'Mcduff', 'Belly', 'Mathew', 'Mark', 'Luke', 'John',
    'Aslam', 'Ham', 'Shem', 'Joshua', 'Jacob', 'Esaw', 'Omar', 'Enoch',
    'Obadia', 'Daniel', 'Samuel', 'Robbo', 'Joebed', 'Ismael', 'Isreal', 'Isabel',
    'Isarat', 'Monk', 'Blip', 'Bacon', 'Danube', 'Friend', 'Darryl', 'Izzy',
    'Crosby', 'Stills', 'Nash', 'Young', 'Cheese', 'Salami', 'Prawn', 'Radish',
    'Egbert', 'Edwy', 'Edgar', 'Edwin', 'Edred', 'Eggpie', 'Bros', 'Sonic',
    'Ziggy', 'Alfred', 'Siggy', 'Hilda', 'Snell', 'Sparks', 'Spook', 'Topcat',
    'Benny', 'Dibble', 'Benker', 'Dosey', 'Beaky', 'Joist', 'Pivot', 'Tree',
    'Bush', 'Grass', 'Seedy', 'Tin', 'Rollo', 'Zippo', 'Nancy', 'Larry',
    'Iggy', 'Nigel', 'Jamie', 'Jesse', 'Leo', 'Virgo', 'Garth', 'Fidel',
    'Idi', 'Che', 'Kirk', 'Spock', 'Maccoy', 'Chekov', 'Uhura', 'Bones',
    'Vulcan', 'Fester', 'Jethro', 'Jimbob', 'Declan', 'Dalek', 'Hickey', 'Chocco',
    'Goch', 'Pablo', 'Renoir', 'Rolf', 'Dali', 'Monet', 'Manet', 'Gaugin',
    'Chagal', 'Kid', 'Hully', 'Robert', 'Piers', 'Raith', 'Jeeves', 'Paster',
    'Adolf', 'Deiter', 'Deni', 'Zark', 'Wizkid', 'Wizard', 'Iain', 'Kitten',
    'Gonner', 'Waster', 'Loser', 'Fodder',
]
