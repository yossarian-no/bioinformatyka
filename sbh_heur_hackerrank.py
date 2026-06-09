#!/usr/bin/env python3
import sys
import random
from collections import defaultdict

sys.setrecursionlimit(200000)

def parse_input(data):
    lines = data.strip().splitlines()
    n = None
    start = None
    k = None
    probes = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('#length'):
            n = int(line.split()[1])
        elif line.startswith('#start'):
            start = line.split()[1].strip().upper()
        elif line.startswith('#probe'):
            k = int(line.split()[1])
        elif not line.startswith('#'):
            probes.append(line.strip().upper())
    return n, start, k, probes

def build_graph(probes):
    graph = defaultdict(list)
    for a in probes:
        suffix = a[1:]
        for b in probes:
            if a != b and b.startswith(suffix):
                graph[a].append(b)
    return graph

def count_covered(dna, probes, k):
    probe_set = set(probes)
    covered = 0
    for i in range(len(dna) - k + 1):
        if dna[i:i+k] in probe_set:
            covered += 1
    return covered

def greedy_build(probes, n, start, graph):
    used = {start}
    dna = start
    current = start
    probe_set = set(probes)

    while len(dna) < n:
        candidates = [nb for nb in graph[current] if nb not in used]
        if candidates:
            # lookahead-1: choose candidate with most unused outgoing neighbors
            best = max(candidates, key=lambda p: sum(
                1 for nb in graph[p] if nb not in used and nb != p
            ))
            dna += best[-1]
            used.add(best)
            current = best
        else:
            # no overlap neighbor found - try any unused probe whose prefix
            # matches current suffix
            suffix = current[1:]
            found = False
            for p in probes:
                if p not in used and p.startswith(suffix):
                    dna += p[-1]
                    used.add(p)
                    current = p
                    found = True
                    break
            if not found:
                dna += 'A'
    return list(dna[:n])

def local_search(dna, probes, n, k, max_iter=1000):
    probe_set = set(probes)
    best = dna[:]
    best_score = count_covered(''.join(best), probes, k)
    bases = ['A', 'T', 'G', 'C']

    for _ in range(max_iter):
        if best_score == len(probes):
            break
        pos = random.randint(0, n - 1)
        old = best[pos]
        for nt in bases:
            if nt == old:
                continue
            best[pos] = nt
            score = count_covered(''.join(best), probes, k)
            if score >= best_score:
                best_score = score
                break
            else:
                best[pos] = old

    return best

def multi_start(probes, n, start, graph, k, restarts=5):
    best_dna = None
    best_score = -1
    probe_list = probes[:]

    for i in range(restarts):
        random.shuffle(probe_list)
        # rebuild graph with shuffled order for variety
        g = defaultdict(list)
        for a in probe_list:
            suffix = a[1:]
            for b in probe_list:
                if a != b and b.startswith(suffix):
                    g[a].append(b)

        s = start if (start and start in set(probe_list)) else probe_list[0]
        dna = greedy_build(probe_list, n, s, g)
        dna = local_search(dna, probes, n, k, max_iter=500)
        score = count_covered(''.join(dna), probes, k)

        if score > best_score:
            best_score = score
            best_dna = dna[:]

        if best_score == len(probes):
            break

    return best_dna or list('A' * n)

def solve(data):
    random.seed(42)
    n, start, k, probes = parse_input(data)

    if not probes:
        return 'A' * (n or 1)
    if k is None:
        k = len(probes[0])
    if n is None:
        n = len(probes) + k - 1

    graph = build_graph(probes)
    s = start if (start and start in set(probes)) else probes[0]

    result = multi_start(probes, n, s, graph, k, restarts=5)
    return ''.join(result[:n])

def main():
    data = sys.stdin.read()
    print(solve(data))

if __name__ == '__main__':
    main()
