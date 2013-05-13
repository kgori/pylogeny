#!/usr/bin/env python

newick = \
    '((a:12, (b:8, c:11):9):10, ((d:16, e:17):5, f:13):6, (g:15, h:14):7):18;'

tree = dpy.Tree.get_from_string(newick, 'newick')


def get_inner_edges(tree):
    inner_edges = [e for e in tree.preorder_edge_iter() if e.is_internal()
                   and e.head_node and e.tail_node]
    return inner_edges


def get_children(inner_edge):
    """ N1: Dendropy trees are built from seed node outwards. If the inner edge
    has the seed node as its tail node, then its children will include the head
    node, which is not what we want.  """
    head_node = inner_edge.head_node
    tail_node = inner_edge.tail_node
    head_children = set(head_node.child_nodes())
    tail_children = set(tail_node.child_nodes()) - set([head_node])  # See N1
    if tail_node.parent_node:
        tail_children |= set([tail_node.parent_node])
    return {'head': head_children, 'tail': tail_children}


def switch(head_node, tail_node):
    pass
