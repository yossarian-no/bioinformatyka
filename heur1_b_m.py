#!/usr/bin/env python3
"""
heur1_b_m.py — Эвристический алгоритм SBH: Binary Chip, Negative Errors
Использование: python heur1_b_m.py < input.xml
"""
import sys, random
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
BASES = ['A','T','G','C']

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

def probe_covered(dna, pos, probe):
    return pos + len(probe) <= len(dna) and all(match(dna[pos+i], ch) for i, ch in enumerate(probe))

def count_covered(dna, probes):
    return sum(1 for p in probes if any(probe_covered(dna, i, p) for i in range(len(dna)-len(p)+1)))

def greedy_build(probes, n, k):
    if not probes:
        return list('ACGT' * (n//4+1))[:n]
    starts = set(probes)
    for p in probes:
        for q in probes:
            if p != q and overlap(q, p):
                starts.discard(p); break
    start = next(iter(starts), probes[0])
    dna = list(expand(start)[0]) if expand(start) else [random.choice(BASES)]*k
    used, current = {start}, start
    while len(dna) < n:
        cands = [p for p in probes if p not in used and overlap(current, p)]
        if cands:
            best = max(cands, key=lambda p: sum(1 for q in probes if q not in used and q != p and overlap(p, q)))
            dna.append(new_nt(best)[0]); used.add(best); current = best
        else:
            rem = [p for p in probes if p not in used]
            if rem:
                suf = ''.join(dna[-(k-1):])
                ns = next((p for p in rem if all(match(suf[i], p[i]) for i in range(min(len(suf),len(p))))), None)
                if ns:
                    used.add(ns); dna.append(new_nt(ns)[0]); current = ns
                else:
                    dna.append(random.choice(BASES))
            else:
                dna.append(random.choice(BASES))
    return dna[:n]

def local_search(dna, probes, n, max_iter=1000):
    best, score = dna[:], count_covered(dna, probes)
    for _ in range(max_iter):
        if score == len(probes): break
        pos = random.randint(0, n-1)
        old = best[pos]
        bs = BASES[:]; random.shuffle(bs)
        for nt in bs:
            if nt == old: continue
            best[pos] = nt
            ns = count_covered(best, probes)
            if ns >= score:
                score = ns; break
            else:
                best[pos] = old
    return best

def multi_start(probes, n, k, restarts=5):
    best_dna, best_score = None, -1
    sh = probes[:]
    for _ in range(restarts):
        random.shuffle(sh)
        dna = greedy_build(sh, n, k)
        dna = local_search(dna, probes, n, max_iter=500)
        sc = count_covered(dna, probes)
        if sc > best_score:
            best_score, best_dna = sc, dna[:]
        if best_score == len(probes): break
    return best_dna or list('A'*n)

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
    random.seed(42)
    n, k, probes = parse_xml(sys.stdin.read())
    print(''.join(multi_start(probes, n, k)[:n]))

if __name__ == '__main__':
    main()
