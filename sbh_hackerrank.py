#!/usr/bin/env python3
import sys
from collections import defaultdict

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

def reconstruct(probes, n, start, graph):
    def dfs(node, dna, visited):
        if len(dna) == n:
            return dna
        if len(dna) > n:
            return None
        for nb in graph[node]:
            if nb in visited:
                continue
            visited.add(nb)
            result = dfs(nb, dna + nb[-1], visited)
            if result is not None:
                return result
            visited.discard(nb)
        return None

    visited = {start}
    return dfs(start, start, visited)

def greedy_fallback(probes, n, start):
    used = {start}
    dna = start
    current = start
    while len(dna) < n:
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
    return dna[:n]

def solve(data):
    n, start, k, probes = parse_input(data)

    if not probes:
        return 'A' * (n or 1)
    if k is None:
        k = len(probes[0])
    if n is None:
        n = len(probes) + k - 1

    graph = build_graph(probes)

    if start and start in set(probes):
        result = reconstruct(probes, n, start, graph)
        if result:
            return result
        return greedy_fallback(probes, n, start)
    else:
        for s in probes:
            result = reconstruct(probes, n, s, graph)
            if result:
                return result
        return greedy_fallback(probes, n, probes[0])

def main():
    data = sys.stdin.read()
    print(solve(data))

if __name__ == '__main__':
    main()
