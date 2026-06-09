#!/usr/bin/env python3

import sys
import xml.etree.ElementTree as ET
from itertools import product

COMPAT = {
    'A': {'A'}, 'T': {'T'}, 'G': {'G'}, 'C': {'C'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'},
    'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'T', 'G', 'C'},
    'Z': {'A', 'T', 'G', 'C'},
    'P': {'A', 'T', 'G', 'C'},
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

def find_start_probe(probes, start_seq):
    # Find probe whose suffix (k-1) is compatible with prefix (k-1) of start_seq
    k = len(probes[0]) if probes else 0
    prefix = start_seq[:k-1]
    for p in probes:
        suffix = p[1:]
        if len(suffix) == len(prefix) and all(match(suffix[i], prefix[i]) for i in range(len(prefix))):
            return p
    return None

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

def solve(probes, n, k, start_seq=None):
    if not probes:
        return 'A' * n
    graph = build_graph(probes)

    # Build ordered list: start probe first, then the rest
    start_probe = None
    if start_seq:
        start_probe = find_start_probe(probes, start_seq)

    ordered = []
    if start_probe:
        ordered.append(start_probe)
    ordered += [p for p in probes if p != start_probe]

    for start in ordered:
        res = backtrack(graph, probes, n, start)
        if res is not None:
            return ''.join(res)
    return ''.join(greedy(probes, n, k))

def parse_xml(data):
    root = ET.fromstring(data)
    n = int(root.get('length', root.get('n', 200)))
    k = int(root.get('k', 10))
    start_seq = root.get('start', None)
    probes = []
    for tag in ('cell', 'probe'):
        for el in root.iter(tag):
            v = el.get('pattern') or el.get('value') or el.text
            if v:
                probes.append(v.strip().upper())
    return n, k, start_seq, list(set(probes))

def main():
    import os
    xml_file = "input.xml"
    if os.path.exists(xml_file):
        with open(xml_file, encoding="utf-8") as f:
            data = f.read()
    else:
        data = sys.stdin.read()

    n, k, start_seq, probes = parse_xml(data)
    result = solve(probes, n, k, start_seq)[:n]
    print(result)

    with open("output_exact.txt", "w") as out:
        out.write(result + "\n")
    print(f"[INFO] Result saved to output_exact.txt (length: {len(result)})")

if __name__ == '__main__':
    main()
