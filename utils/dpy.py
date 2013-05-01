#!/usr/bin/env python

import dendropy


# DENDROPY UTILS

def bifurcate_base(newick):
    t = dendropy.Tree.get_from_string(newick, 'newick')
    t.resolve_polytomies()
    newick_string = t.as_newick_string() + ';\n'
    return newick_string


def trifurcate_base(newick):
    t = dendropy.Tree.get_from_string(newick, 'newick')
    t.deroot()
    return t.as_newick_string() + ';\n'


def brlen_sum(dpy_tree):
    tot = 0
    for n in dpy_tree.preorder_node_iter():
        if n.edge_length:
            tot += n.edge_length
    return tot


def convert_dendropy_to_newick(tree):
    newick = tree.as_newick_string()
    return (newick if newick.endswith(';') else newick + ';')


def convert_to_dendropy_tree(tree, taxon_set=None):
    if taxon_set is None:
        return dendropy.Tree.get_from_string(tree.newick, 'newick')
    return dendropy.Tree.get_from_string(tree.newick, 'newick',
            taxon_set=taxon_set)


def convert_to_dendropy_trees(trees):
    taxa = dendropy.TaxonSet()
    return [convert_to_dendropy_tree(tree, taxa) for tree in trees]


def check_diff_top(check_tree, tree_list):
    """ Returns True if topologies are different, False if they are the same
    (unweighted RF distance = 0) """

    checklist = []
    t1 = convert_to_dendropy_tree(check_tree)
    tree_list = convert_to_dendropy_trees(tree_list)

    for tree in tree_list:
        if tree.symmetric_difference(t1) == 0:
            return False
    return True


def check_rooted(newick):
    if newick is None or newick == '':
        return None
    t = dendropy.Tree.get_from_string(newick, 'newick')
    root_degree = len(t.seed_node.child_nodes())
    return root_degree == 2


def get_rf_distance(dpy_tree1, dpy_tree2):
    return dpy_tree1.symmetric_difference(dpy_tree2)


def get_wrf_distance(dpy_tree1, dpy_tree2):
    return dpy_tree1.robinson_foulds_distance(dpy_tree2)


def get_euc_distance(dpy_tree1, dpy_tree2):
    return dpy_tree1.euclidean_distance(dpy_tree2)

def ntaxa(tree):
    t = convert_to_dendropy_tree(tree)
    return len(t.leaf_nodes())

def taxon_labels(dpy_tree):
    return set([n.taxon.label for n in dpy_tree.leaf_nodes()])

def length(tree):
    t = convert_to_dendropy_tree(tree)
    return t.length()

def pruned_pair(tree1, tree2):
    t1 = convert_to_dendropy_tree(tree1)
    t2 = convert_to_dendropy_tree(tree2)
    tax1 = taxa(t1)
    tax2 = taxa(t2)
    in1not2 = list(tax1 - tax2)
    in2not1 = list(tax2 - tax1)
    t1.prune_taxa_with_labels(in1not2)
    t2.prune_taxa_with_labels(in2not1)
    n1 = convert_dendropy_to_newick(t1)
    n2 = convert_dendropy_to_newick(t2)
    p1 = tree1.copy()
    p2 = tree2.copy()
    p1.newick = n1
    p2.newick = n2
    return p1, p2

def print_plot(tree):
    return convert_to_dendropy_tree(tree).print_plot()

