#!/usr/bin/env python

class TreeEdgeError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg

def rootcheck(edge, msg='This is the root edge'):
    if not edge.tail_node:
        raise TreeEdgeError(msg)

def edge_length_check(length, edge):
    try:
        assert 0 <= length <= edge.length
    except AssertionError:
        if length < 0:
            raise TreeEdgeError('Negative edge-lengths are disallowed')
        raise TreeEdgeError(
            'This edge isn\'t long enough to prune at length {0}\n'
            '(Edge length = {1})'.format(length, edge.length))    

def prune(tree, edge, length=None):

    length = length or edge.length
    edge_length_check(length, edge)

    n = edge.head_node
    t.prune_subtree(n)
    n.edge_length = length
    return n

def regraft(tree, edge, node, length=None):
    length = length or edge.length/2. # Length measured from head to tail
    edge_length_check(length, edge)
    rootcheck(edge, 'SPR regraft is not allowed on the root edge')

    t = edge.tail_node
    h = edge.head_node
    new = t.new_child(edge_length=edge.length-length)
    t.remove_child(h)
    new.add_child(h, edge_length=length)
    new.add_child(node)
    tree.update_splits()

def are_sibling_edges(e1, e2):
    e1check = e1.head_node in e2.head_node.sister_nodes()
    e2check = e2.head_node in e1.head_node.sister_nodes()
    if e1check != e2check:
        raise TreeEdgeError(
            'Edges {0} and {1} have an inconsistent sibling relationship')
    return e1check and e2check

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



def rspr(tree, inplace=False, disallow_sibling_sprs=False):
    
    if not inplace:
        tree = tree.copy()
    else:
        tree = tree
    
    tree.print_plot()
    edge_names, node_names, r_edges, r_nodes = _name_things(tree)

    edges = [e for e in tree.preorder_edge_iter()
             if e.head_node and e.tail_node]
    pruning_edge = random.choice(edges)
    
    if disallow_sibling_sprs:
        for e in pruning_edge.adjacent_edges:
            edges.remove(e)

    for n in pruning_edge.head_node.preorder_iter():
        edges.remove(n.edge)

    regrafting_edge = random.choice(edges)
    length1 = random.uniform(0, pruning_edge.length)
    length2 = random.uniform(0, regrafting_edge.length)

    print ('Pruning edge = {0},{1}\n'
        'Regrafting edge = {2},{3}'.format(edge_names[pruning_edge],
            length1, 
            edge_names[regrafting_edge],
            length2))
    spr(tree, pruning_edge, regrafting_edge, length1, length2)
    tree.print_plot()
    return tree, len(tree)


def spr(tree, pruning_edge, regrafting_edge, length1=None, length2=None):

    try:
        assert (regrafting_edge is not pruning_edge)
    except AssertionError:
        raise TreeEdgeError('Can\'t regraft on the pruning edge')

    try:
        descendents = [n.edge for n in pruning_edge.head_node.preorder_iter()]
        assert regrafting_edge not in descendents
    except AssertionError:
        raise TreeEdgeError('Regraft edge is on the pruned subtree')

    sister_nodes = pruning_edge.head_node.sister_nodes()
    if regraft == pruning_edge.tail_node.edge and len(sister_nodes) == 1:
        length2 += sister_nodes[0].length
    n = prune(tree, pruning_edge, length1)
    

    regraft(tree, regrafting_edge, n, length2)
