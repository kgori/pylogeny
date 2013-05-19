#!/usr/bin/env python

import dendropy as dpy
import random


def get_inner_edges(tree):
    inner_edges = [e for e in tree.preorder_edge_iter() if e.is_internal()
                   and e.head_node and e.tail_node]
    return inner_edges


def get_children(tree, inner_edge):
    """ N1: Edges are directional in dendropy trees. The head node of an edge is
    automatically a child of the tail node, but we don't want this. """

    h = inner_edge.head_node
    t = inner_edge.tail_node
    if not tree.seed_node == t:
        original_seed = tree.seed_node
        tree.reseed_at(t)
    else: 
        original_seed = None
    head_children = h.child_nodes()
    tail_children = list(set(t.child_nodes()) - set([h])) # See N1
    if original_seed:
        tree.reseed_at(original_seed)
    return {'head': head_children, 'tail': tail_children}

def nni(
    tree,
    edge,
    nn_head,
    nn_tail,
    ):
    original_seed = tree.seed_node
    h = edge.head_node
    t = edge.tail_node
    tree.reseed_at(t)
    assert nn_head.parent_node == h
    assert nn_tail.parent_node == t
    h.remove_child(nn_head)
    t.remove_child(nn_tail)
    h.add_child(nn_tail)
    t.add_child(nn_head)
    tree.reseed_at(original_seed, update_splits=True)

def _name_things(tree):
    edges = {}
    nodes = {None: 'root'}
    for n in tree.postorder_node_iter():
        nodes[n] = '.'.join([str(x.taxon) for x in n.leaf_nodes()])
    for e in tree.preorder_edge_iter():
        edges[e] = ' ---> '.join([nodes[e.tail_node], nodes[e.head_node]])

    r_edges = {value: key for key,value in edges.items()}
    r_nodes = {value: key for key,value in nodes.items()}
    return edges, nodes, r_edges, r_nodes

def rnni(tree, verbosity=1):
    if verbosity > 0:
        # name things
        edge_names, node_names, r_edges, r_nodes = _name_things(tree)
    
    edges = get_inner_edges(tree)
    e = random.choice(edges)
    children = get_children(tree, e)
    (h, t) = (random.choice(children['head']),
              random.choice(children['tail']))
    if verbosity > 0:
        print 'Selected edge [ {0} ]'.format(edge_names[e]) 
        print 'Interchanging neighbours {0} and {1}'.format(node_names[h], 
            node_names[t])
    nni(tree, e, h, t)
    return tree

def rnni_short_version(tree):
    e = random.choice(get_inner_edges(tree))
    children = get_children(tree, e)
    (h, t) = (random.choice(children['head']),
              random.choice(children['tail']))
    nni(tree, e, h, t)
    return tree

if __name__ == '__main__':

    newick = \
        '((a:12, (b:8, c:11):9):10, ((d:16, e:17):5, f:13):6, (g:15, h:14):7):18;'

    tree = dpy.Tree.get_from_string(newick, 'newick')
    original_seed = tree.seed_node

    edges, nodes, r_edges, r_nodes = _name_things(tree)

    e = r_edges['d.e.f ---> d.e']
    nbh = r_nodes['d']
    nbt = tree.seed_node
    tree.print_plot()
    print tree.as_newick_string()
    nni(tree, e, nbh, nbt)
    tree.print_plot()
    print tree.as_newick_string()
