#!/usr/bin/env python

import random
import re
import dendropy as dpy
from numpy.random import gamma
from ..errors import FileError, filecheck
from ..utils import dpy as utils_dpy


def blob(master, ntimes, ntrees):
    """ Placeholder for a sklearn.datasets.make_blobs style function 
    to generate a group of similar trees from a master tree """
    blob = []
    for tree in range(ntrees):
        manip = TreeManip(master)
        for time in range(ntimes):
            manip.spr()
        blob.append(manip.tree)
    return blob

class TreeManip(object):

    def __init__(self, tree):
        if not isinstance(tree, Tree):
            raise TypeError('TreeManip class should be initialised with \'Tree\' object. Got'
                            , type(tree))
        self.tree = tree.copy()


    def __str__(self):
        return 'TreeManip object with tree:\n' + str(self.tree)

    def convert_to_dendropy_tree(self):
        """Takes Tree object, returns dendropy.Tree object"""

        return utils_dpy.convert_to_dendropy_tree(self.tree)

    def dendropy_as_newick(self, dpy_tree):
        return utils_dpy.convert_dendropy_to_newick(dpy_tree)

    def randomise_branch_lengths(
        self,
        i=(1, 1),
        l=(1, 1),
        distribution_func=gamma,
        output_format=5,
        ):
        """ i and l are tuples describing the parameters of the distribution
        function for inner edges and leaves, respectively. distribution_func is
        a function generating samples from a probability distribution (eg gamma,
        normal ...) """

        t = self.convert_to_dendropy_tree()
        for n in t.preorder_node_iter():
            if n.is_internal():
                n.edge.length = max(0, distribution_func(*i))
            else:
                n.edge.length = max(0, distribution_func(*l))
        tree_copy = self.tree.copy()
        tree_copy.newick = self.dendropy_as_newick(t)
        self.tree = tree_copy
        return self.tree

    def prune_to_subset(self, subset):
        t = utils_dpy.convert_to_dendropy_tree(self.tree)
        labels = utils_dpy.taxon_labels(t)
        if not subset.issubset(labels):
            print '"subset" is not a subset'
        else: 
            t.retain_taxa_with_labels(subset)
            self.tree.newick = self.dendropy_as_newick(t)
        return self.tree


    def randomise_labels(self):
        t = self.convert_to_dendropy_tree()
        names = [l.taxon.label for l in t.leaf_iter()]
        names_copy = names[:]
        random.shuffle(names_copy)
        for l in t.leaf_iter():
            l.taxon.label = names_copy.pop()
        tree_copy = self.tree.copy()
        tree_copy.newick = self.dendropy_as_newick(t)
        self.tree = tree_copy
        return self.tree

    def scale(self, scaling_factor):

        t = self.convert_to_dendropy_tree()
        for e in t.preorder_edge_iter():
            if e.length:
                e.length *= scaling_factor
        tree_copy = self.tree.copy()
        tree_copy.newick = self.dendropy_as_newick(t)
        self.tree = tree_copy
        return self.tree

    def strip(self):

        t = self.convert_to_dendropy_tree()
        for e in t.preorder_edge_iter():
            e.length = None
        tree_copy = self.tree.copy()
        tree_copy.newick = self.dendropy_as_newick(t)
        self.tree = tree_copy
        return self.tree

    def nni(self):
        """Function to perform a random nearest-neighbour interchange on a tree
        using Dendropy
        
        The dendropy representation of a tree is as if rooted (even when it's
        not) The node which acts like the root is called the seed node, and this
        can sit in the middle of an edge which would be a leaf edge in the
        unrooted tree. NNI moves are only eligible on internal edges, so we need
        to check if the seed node is sat on a real internal edge, or a fake
        one."""

        tree = self.convert_to_dendropy_tree()
        seed = tree.seed_node
        resolved = False                    
        if len(seed.child_nodes()) > 2:     # If root is multifurcating we
            print 'Resolve root trisomy'    # need to force it to bifurcate,
            tree.resolve_polytomies()       # and re-establish multifurcation
            resolved = true                 # later

        # Make a list of internal edges not including the root edge
        edge_list = list(tree.preorder_edge_iter(lambda edge: \
                         (True if edge.is_internal() and edge.head_node != seed
                         and edge.tail_node != seed else False)))

        # Test whether root edge is eligible (i.e., the edge is not a
        # leaf when the tree is unrooted). If the test passes, add 'root'
        # to the list of eligible edges
        if not any([x.is_leaf() for x in seed.child_nodes()]):
            edge_list += ['root']

        chosen_edge = random.choice(edge_list)  # Choose the edge around which
        print chosen_edge                       # to do the NNI

        # Reroot at chosen edge, if necessary
        if chosen_edge != 'root':
            tree.reroot_at_edge(chosen_edge, length1=chosen_edge.length / 2,
                                length2=chosen_edge.length / 2,
                                delete_outdegree_one=False)
            root = tree.seed_node
        else:
            root = seed

        # To do the swap: find the nodes on either side of root
        (child_left, child_right) = root.child_nodes()

        # Pick a child node from each of these
        neighbour1 = random.choice(child_left.child_nodes())
        neighbour2 = random.choice(child_right.child_nodes())

        # Prune the chosen nearest neighbours - but don't
        # do any tree structure updates
        tree.prune_subtree(neighbour1, update_splits=False,
                           delete_outdegree_one=False)
        tree.prune_subtree(neighbour2, update_splits=False,
                           delete_outdegree_one=False)

        # Reattach the pruned neighbours to the opposite side
        # of the tree
        child_left.add_child(neighbour2)
        child_right.add_child(neighbour1)

        # Reroot the tree using the original seed node, and
        # update splits
        if not chosen_edge == 'root':
            tree.reroot_at_node(seed, update_splits=True)
        else:
            tree.update_splits()

        if resolved:
            print 'Reinstating root trisomy'
            tree.deroot()

        newick = self.dendropy_as_newick(tree)
        if tree.is_rooted:
            newick = '[&R] ' + newick

        tree_copy = self.tree.copy() 
        tree_copy.newick = newick
        self.tree = tree_copy
        return self.tree

    

class Tree(object):

    """ Class for storing the results of phylogenetic inference """

    name_regex = re.compile('([A-Za-z0-9\-_]+).([A-Za-z0-9\-_]+)(?=_phyml_)')
    score_regex = re.compile('(?<=Log-likelihood: ).+')

    def __init__(
        self,
        newick=None,
        score=0,
        output=None,
        program=None,
        name=None,
        rooted=None,
        ):

        self.newick = newick
        self.score = score
        self.output = output
        self.program = program
        self.name = name
        self.rooted = utils_dpy.check_rooted(newick)

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


    def __eq__(self, other):

        return all([self.name == other.name, self.newick == other.newick,
                   self.program == other.program, self.score == other.score,
                   self.output == other.output])


    def __len__(self):
        return utils_dpy.ntaxa(self)

    def __and__(self, other):
        return self.intersection(other)

    def copy(self):
        copy = self.__new__(type(self))
        copy.__dict__ = {key: value for (key, value) in self.__dict__.items()}
        return copy


    def labels(self):
        t = utils_dpy.convert_to_dendropy_tree(self)
        return utils_dpy.taxon_labels(t)

    def intersection(self, other):
        taxa1 = self.labels()
        taxa2 = other.labels()
        return taxa1 & taxa2

    def pruned_pair(self, other):
        common = self & other
        p1 = self.prune_to_subset(common)
        p2 = other.prune_to_subset(common)
        return (p1, p2)

    def prune_to_subset(self, subset):
        t = TreeManip(self)
        tree = t.prune_to_subset(subset)
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

        writer = open(outfile, 'w')
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
        writer.close()
        return outfile

    def load_phyml_strings(
        self,
        tree,
        stats,
        name=None,
        program='phyml',
        ):

        score = float(self.score_regex.search(stats).group(1))
        self.program = program
        self.newick = tree
        self.output = stats
        self.score = score
        self.name = name
        self.rooted = utils_dpy.check_rooted(tree)

    def load_phyml_files(
        self,
        tree_file,
        stats_file,
        name=None,
        program='phyml',
        ):
        """ Loads phyml results into existing tree object - returns None """

        exit = False
        for f in (tree_file, stats_file):
            try:
                filecheck_and_raise(f)
            except FileError, e:
                print e
                exit = True

        if exit:
            print 'Results were not loaded'
            raise FileError()

        if not name:
            name = self.name_regex.search(tree_file).group(1)
        newick = open(tree_file).read()
        stats = open(stats_file).read()
        self.load_phyml_strings(newick, stats, name=name, program=program)

    @classmethod
    def new_tree_from_phyml_results(
        cls,
        tree_file,
        stats_file,
        program='phyml',
        ):
        """ classmethod version of load_phyml_files - returns a new Tree object
        """

        new_tree = cls()
        new_tree.load_phyml_files(tree_file, stats_file, program=program)
        return new_tree

    @classmethod
    def new_tree_from_phyml_strings(
        cls,
        tree,
        stats,
        program='phyml',
        ):

        new_tree = cls()
        new_tree.load_phyml_strings(tree, stats, program=program)
        return new_tree

    def scale(self, scale_factor):
        return TreeManip(self).scale(scale_factor)

    def strip(self):
        return TreeManip(self).strip()

    def length(self):
        return utils_dpy.length(self)

    def ntaxa(self):
        return utils_dpy.ntaxa(self)

    def print_plot(self):
        utils_dpy.print_plot(self)

    def rfdist(self, other):
        s = utils_dpy.convert_to_dendropy_tree(self)
        o = utils_dpy.convert_to_dendropy_tree(other)
        return utils_dpy.get_rf_distance(s, o)

    def wrfdist(self, other):
        s = utils_dpy.convert_to_dendropy_tree(self)
        o = utils_dpy.convert_to_dendropy_tree(other)
        return utils_dpy.get_wrf_distance(s, o)

    def eucdist(self, other):
        s = utils_dpy.convert_to_dendropy_tree(self)
        o = utils_dpy.convert_to_dendropy_tree(other)
        return utils_dpy.get_euc_distance(s, o)

    @classmethod
    def new_yule(
        self,
        nspecies,
        names=None,
        cf=False,
        ):
        g = TreeGen(nspecies, names, cf=cf)
        return g.yule()

    @classmethod
    def new_coal(
        self,
        nspecies,
        names=None,
        cf=False,
        ):
        g = TreeGen(nspecies, names, cf=cf)
        return g.coal()

    @classmethod
    def new_rtree(
        self,
        nspecies,
        names=None,
        cf=False,
        ):
        g = TreeGen(nspecies, names, cf=cf)
        return g.rtree()

    def gene_tree(self, scale_to=None):
        """ Returns a constrained Kingman gene tree using self as species 
        tree. Optionally rescales to height given by scale_to parameter """
        g = TreeGen(template=self)
        return g.gene_tree(scale_to)['gene_tree']


class TreeGen(object):

    def __init__(
        self,
        nspecies=16,
        names=None,
        template=None,
        cf=False,
        ):

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
        taxon_set = dpy.TaxonSet(self.names)
        tree = dpy.treesim.pure_kingman(taxon_set)
        newick = '[&R] ' + utils_dpy.convert_dendropy_to_newick(tree)
        return Tree(newick)

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
        tree = utils_dpy.convert_to_dendropy_tree(t)

        for leaf in tree.leaf_iter():
            leaf.num_genes = 1

        tree_height = tree.seed_node.distance_from_root() \
            + tree.seed_node.distance_from_tip()

        if scale_to:
            population_size = tree_height / scale_to

        for edge in tree.preorder_edge_iter():
            edge.pop_size = population_size

        gene_tree = dpy.treesim.constrained_kingman(tree)[0]

        if trim_names:
            for leaf in gene_tree.leaf_iter():
                leaf.taxon.label = leaf.taxon.label.replace('\'', '').split('_'
                        )[0]

        newick = '[&R] ' + utils_dpy.convert_dendropy_to_newick(gene_tree)

        return {'gene_tree': Tree(newick), 'species_tree': t}

    def rtree(self):
        m = TreeManip(self.yule())
        m.randomise_labels()
        return m.randomise_branch_lengths()

    def yule(self):
        taxon_set = dpy.TaxonSet(self.names)
        tree = dpy.treesim.uniform_pure_birth(taxon_set)
        newick = '[&R] ' + utils_dpy.convert_dendropy_to_newick(tree)
        return Tree(newick)


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
