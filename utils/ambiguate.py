#!/usr/bin/env python
import treeCl
from treeCl.lib.remote.datastructs.seq import Seq

docstring = '''
Combine two unphased sequences into a single sequence with ambiguity codes
Usage: ambiguate.py <FILENAME=eg:phylip.phy>
Output: <phylip_ambig.phy>

'''

ambiguities = {
    'a': frozenset(['a']),
    'c': frozenset(['c']),
    'g': frozenset(['g']),
    't': frozenset(['t']),
    'y': frozenset(['c', 't']),
    'r': frozenset(['a', 'g']),
    'w': frozenset(['a', 't']),
    's': frozenset(['g', 'c']),
    'k': frozenset(['t', 'g']),
    'm': frozenset(['c', 'a']),
    'd': frozenset(['a', 'g', 't']),
    'v': frozenset(['a', 'c', 'g']),
    'h': frozenset(['a', 'c', 't']),
    'b': frozenset(['c', 'g', 't']),
    'x': frozenset(['a', 'c', 'g', 't']),
    'n': frozenset(['a', 'c', 'g', 't']),
}

ambiguities_rev = {v:k for (k,v) in ambiguities.items()}
ambiguities_rev[frozenset(['a', 'c', 'g', 't'])] = 'n'

def get_prefixes(r):
    prefixes = list()
    for x in r.headers:
        prefix = '.'.join(x.split('.')[:-1])
        if not prefix in prefixes:
            prefixes.append(prefix)
    return prefixes

def get_ambiguity(a, b):
    upper = False
    if a.isupper():
        upper = True
        a = a.lower()
    if b.isupper():
        upper = True
        b = b.lower()
    s1 = ambiguities[a]
    s2 = ambiguities[b]
    union = s1 | s2
    ambig = ambiguities_rev[union]
    
    return ambig.upper() if upper else ambig

def ambiguate(seq1, seq2):
    combination = list()
    z = zip(seq1, seq2)
    for (a,b) in z:
        if a == b:
            combination.append(a)
        else:
            if a == '-' or b == '-':
                combination.append('-')
            else:
                ambig = get_ambiguity(a, b)
                combination.append(ambig)
    return ''.join(combination)

def get_seqs(rec, pref):
    return rec.mapping[pref+'.1'], rec.mapping[pref+'.2']

if __name__ == '__main__':
    import sys
    if len(sys.argv) == 1:
        print 'No file entered - syntax = ambiguate.py <FILE>'
        sys.exit(1)
    f = sys.argv[1]
    if f == '-h' or f == '--help':
        print docstring
        sys.exit()
    rec = Seq(f, 'phylip')
    prefixes = get_prefixes(rec)
    print prefixes
    headers = list()
    sequences = list()
    for pref in prefixes:
        print pref
        seq1, seq2 = get_seqs(rec, pref)
        combin = ambiguate(seq1, seq2)
        headers.append(pref)
        sequences.append(combin)
    newrec = Seq(headers=headers, sequences=sequences)
    out = f[:f.rindex('.phy')] + '_ambig.phy'
    newrec.write_phylip(out, interleaved=True)
