#!/usr/bin/env python3
"""
exact_b_m.py — Точный алгоритм SBH: Binary Chip, Negative Errors
Использование: python exact_b_m.py < input.xml
"""
import sys
import xml.etree.ElementTree as ET
from itertools import product

COMPAT = {
    'A': {'A'}, 'T': {'T'}, 'G': {'G'}, 'C': {'C'},
    'R': {'A','G'}, 'Y': {'C','T'}, 'S': {'G','C'}, 'W': {'A','T'},
    'K': {'G','T'}, 'M': {'A','C'},
    'B': {'C','G','T'}, 'D': {'A','G','T'},
    'H': {'A','C','T'}, 'V': {'A','C','G'},
    'N': {'A','T','G','C'},
}

def match(a, b):
    return bool(COMPAT.get(a, {a}) & COMPAT.get(b, {b}))

def overlap(pa, pb, shift=1):
    ov = len(pa) - shift
    return ov <= 0 or all(match(pa[-ov+i], pb[i]) for i in range(ov))

def expand(probe):
    opts = [sorted(COMPAT.get(c, {c})) for c in probe]
    return [''.join(x) for x in product(*opts)]

def new_nt(probe):
    return sorted(COMPAT.get(probe[-1], {probe[-1]}))

def build_graph(probes):
    g = {p: [] for p in probes}
    for a in probes:
        for b in probes:
            if a != b and overlap(a, b):
                g[a].append(b)
    return g

def backtrack(graph, probes, n, start):
    for start_seq in expand(start):
        path, dna, vis = [start], list(start_seq), {start}

        def dfs():
            if len(dna) == n:
                return True
            if len(dna) > n:
                return False
            for nb in graph.get(path[-1], []):
                if nb in vis:
                    continue
                for nt in new_nt(nb):
                    path.append(nb); dna.append(nt); vis.add(nb)
                    if dfs():
                        return True
                    path.pop(); dna.pop(); vis.discard(nb)
            return False

        if dfs():
            return dna
    return None

def greedy(probes, n, k):
    if not probes:
        return list('A' * n)
    current = probes[0]
    dna = list(expand(current)[0]) if expand(current) else ['A'] * k
    used = {current}
    while len(dna) < n:
        best = next((p for p in probes if p not in used and overlap(current, p)), None)
        if best:
            dna.append(new_nt(best)[0]); used.add(best); current = best
        else:
            dna.append('A')
    return dna[:n]

def solve(probes, n, k):
    if not probes:
        return 'A' * n
    graph = build_graph(probes)
    for start in probes:
        res = backtrack(graph, probes, n, start)
        if res is not None:
            return ''.join(res)
    return ''.join(greedy(probes, n, k))

def parse_xml(data):
    root = ET.fromstring(data)
    n = int(root.get('length', root.get('n', 200)))
    k = int(root.get('k', 10))
    probes = []
    for tag in ('cell', 'probe'):
        for el in root.iter(tag):
            v = el.get('pattern') or el.get('value') or el.text
            if v:
                probes.append(v.strip().upper())
    return n, k, list(set(probes))

def main():
    n, k, probes = parse_xml(sys.stdin.read())
    print(solve(probes, n, k)[:n])

if __name__ == '__main__':
    main()
