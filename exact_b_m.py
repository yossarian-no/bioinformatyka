#!/usr/bin/env python3

import sys
import xml.etree.ElementTree as ET

try:
    import pulp
except ImportError:
    sys.stdout.write("\n")
    sys.exit(0)


COMPAT = {
    'A': {'A'}, 'T': {'T'}, 'G': {'G'}, 'C': {'C'},

    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},

    'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'V': {'A', 'C', 'G'},

    'N': {'A', 'C', 'G', 'T'},
    'Z': {'A', 'C', 'G', 'T'},
    'P': {'A', 'C', 'G', 'T'},
}

BASE_ORDER = ['A', 'C', 'G', 'T']
START = -1

# Для сдачи с лимитом 100 секунд лучше оставить 95.
# Для настоящего exact без ограничения времени поставь None.
TIME_LIMIT_SECONDS = 95


def match(a, b):
    return bool(COMPAT.get(a, {a}) & COMPAT.get(b, {b}))


def concretize_symbol(ch):
    options = COMPAT.get(ch, {ch})

    for base in BASE_ORDER:
        if base in options:
            return base

    return 'A'


def concretize_sequence(seq):
    return ''.join(concretize_symbol(ch) for ch in seq)


def parse_xml(data):
    root = ET.fromstring(data)

    n = int(root.get('length', root.get('n', 0)))

    start_seq = root.get('start', '')
    start_seq = start_seq.strip().upper() if start_seq else ''

    probes = []

    for el in root.iter('cell'):
        if el.text and el.text.strip():
            probes.append(el.text.strip().upper())

    probes = sorted(set(probes))

    if probes:
        k = len(probes[0])
    else:
        k = int(root.get('k', 10))

    if n <= 0:
        if start_seq:
            n = len(start_seq)
        else:
            n = k

    return n, k, start_seq, probes


def get_min_shift(u, v):
    """
    cost(u, v) = k - overlap(u, v)

    Возвращает количество новых нуклеотидов, которое нужно добавить,
    если поставить v после u.
    """
    k = len(v)
    max_overlap = min(len(u), k - 1)

    for shift in range(1, k + 1):
        overlap_len = k - shift

        if overlap_len > max_overlap:
            continue

        if overlap_len <= 0:
            return k

        start = len(u) - overlap_len
        ok = True

        for i in range(overlap_len):
            if not match(u[start + i], v[i]):
                ok = False
                break

        if ok:
            return shift

    return k


def append_probe(dna, probe, cost, n):
    remaining = n - len(dna)

    if remaining <= 0:
        return dna

    cost = min(cost, remaining)

    if cost <= 0:
        return dna

    addition = ''.join(concretize_symbol(ch) for ch in probe[-cost:])

    return dna + addition


def fallback_sequence(n, start_seq):
    """
    Запасной вывод, чтобы программа всегда завершалась и печатала ДНК.
    ВАЖНО: это не exact-решение. Это только технический fallback.
    """
    if start_seq:
        dna = concretize_sequence(start_seq)
    else:
        dna = ''

    if len(dna) < n:
        dna += 'A' * (n - len(dna))

    return dna[:n]


def build_graph(n, start_seq, probes):
    """
    Строим полный ориентированный граф перекрытий.

    Вершины:
        0..m-1 — зонды
        START = -1 — виртуальный старт

    Дуга i -> j имеет стоимость c_ij.
    """
    m = len(probes)
    arcs = []
    costs = {}

    budget = n - len(start_seq)

    if budget < 0:
        budget = 0

    for j in range(m):
        c = get_min_shift(start_seq, probes[j]) if start_seq else len(probes[j])

        if c <= budget:
            arcs.append((START, j))
            costs[(START, j)] = c

    for i in range(m):
        u = probes[i]

        for j in range(m):
            if i == j:
                continue

            v = probes[j]
            c = get_min_shift(u, v)

            if c <= budget:
                arcs.append((i, j))
                costs[(i, j)] = c

    return arcs, costs, budget


def selected_value(var):
    val = pulp.value(var)
    return val is not None and val > 0.5


def find_subtours(x, y, m):
    """
    Находит отдельные циклы в текущем ILP-решении.

    В корректном решении должна быть одна ścieżka wychodząca ze startu.
    Нельзя иметь отдельные циклы вида:
        a -> b -> c -> a
    """
    selected_nodes = set()

    for i in range(m):
        if selected_value(y[i]):
            selected_nodes.add(i)

    successor = {}

    for (i, j), var in x.items():
        if selected_value(var):
            successor[i] = j

    reachable = set()
    cur = successor.get(START)

    while cur is not None and cur not in reachable:
        reachable.add(cur)
        cur = successor.get(cur)

    unreachable = selected_nodes - reachable

    cycles = []
    seen_global = set()

    for node in list(unreachable):
        if node in seen_global:
            continue

        path = []
        pos = {}
        cur = node

        while cur is not None and cur not in pos and cur not in seen_global:
            pos[cur] = len(path)
            path.append(cur)
            seen_global.add(cur)
            cur = successor.get(cur)

        if cur in pos:
            cycle = path[pos[cur]:]

            if len(cycle) > 0:
                cycles.append(cycle)

    return cycles


def build_dna_from_solution(x, costs, probes, start_seq, n):
    dna = concretize_sequence(start_seq)
    current = START
    visited = set()

    while current in x_successor(x):
        successor = x_successor(x)[current]

        if successor in visited:
            break

        visited.add(successor)

        cost = costs[(current, successor)]
        dna = append_probe(dna, probes[successor], cost, n)

        current = successor

        if len(dna) >= n:
            break

    if len(dna) < n:
        dna += 'A' * (n - len(dna))

    return dna[:n]


def x_successor(x):
    succ = {}

    for (i, j), var in x.items():
        if selected_value(var):
            succ[i] = j

    return succ


def solve_exact_ilp(n, k, start_seq, probes):
    """
    Dokładny algorytm oparty o ILP:

    x_ij = 1, jeżeli łuk i -> j należy do rozwiązania
    y_i  = 1, jeżeli sonda i została użyta

    Cel:
        max sum y_i

    Ograniczenia:
        - start ma co najwyżej jeden łuk wychodzący
        - użyty wierzchołek ma dokładnie jeden łuk wchodzący
        - użyty wierzchołek ma co najwyżej jeden łuk wychodzący
        - suma kosztów łuków nie przekracza budżetu długości
        - brak osobnych cykli
    """
    m = len(probes)

    if m == 0:
        return fallback_sequence(n, start_seq)

    if len(start_seq) >= n:
        return concretize_sequence(start_seq)[:n]

    arcs, costs, budget = build_graph(n, start_seq, probes)

    if not arcs:
        return fallback_sequence(n, start_seq)

    model = pulp.LpProblem("Exact_Binary_Chip_Err_Minus_SBH", pulp.LpMaximize)

    x = {
        (i, j): pulp.LpVariable(f"x_{i}_{j}", lowBound=0, upBound=1, cat="Binary")
        for (i, j) in arcs
    }

    y = {
        i: pulp.LpVariable(f"y_{i}", lowBound=0, upBound=1, cat="Binary")
        for i in range(m)
    }

    model += pulp.lpSum(y[i] for i in range(m))

    model += (
        pulp.lpSum(x[(START, j)] for j in range(m) if (START, j) in x) <= 1
    )

    for a in range(m):
        incoming = pulp.lpSum(
            x[(i, a)]
            for i in [START] + list(range(m))
            if (i, a) in x
        )

        outgoing = pulp.lpSum(
            x[(a, j)]
            for j in range(m)
            if (a, j) in x
        )

        model += incoming == y[a]
        model += outgoing <= y[a]

    model += (
        pulp.lpSum(costs[(i, j)] * x[(i, j)] for (i, j) in x) <= budget
    )

    while True:
        if TIME_LIMIT_SECONDS is None:
            solver = pulp.PULP_CBC_CMD(msg=False)
        else:
            solver = pulp.PULP_CBC_CMD(msg=False, timeLimit=TIME_LIMIT_SECONDS)

        status = model.solve(solver)
        status_name = pulp.LpStatus[status]

        if status_name != "Optimal":
            return fallback_sequence(n, start_seq)

        cycles = find_subtours(x, y, m)

        if not cycles:
            break

        for cycle in cycles:
            if len(cycle) == 1:
                node = cycle[0]

                if (node, node) in x:
                    model += x[(node, node)] <= 0

                continue

            cycle_vars = []

            cycle_set = set(cycle)

            for u in cycle:
                for v in cycle:
                    if u != v and (u, v) in x:
                        cycle_vars.append(x[(u, v)])

            if cycle_vars:
                model += pulp.lpSum(cycle_vars) <= len(cycle_set) - 1

    return build_dna_from_solution(x, costs, probes, start_seq, n)


def main():
    data = sys.stdin.read()

    if not data.strip():
        return

    n, k, start_seq, probes = parse_xml(data)

    result = solve_exact_ilp(n, k, start_seq, probes)

    sys.stdout.write(result + "\n")


if __name__ == "__main__":
    main()
