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
}

BASE_ORDER = ['A', 'C', 'G', 'T']

# A limit slightly shorter than the platform limit leaves enough time
# to parse the input and print the result.
TIME_LIMIT_SECONDS = 95
DEBUG = False


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


def fallback_sequence(n, start_seq):
    """
    Return a valid-length DNA sequence when the solver has no solution.

    This is only a technical fallback and is not guaranteed to be optimal.
    """
    if start_seq:
        dna = concretize_sequence(start_seq)
    else:
        dna = ''

    if len(dna) < n:
        dna += 'A' * (n - len(dna))

    return dna[:n]


def selected_value(var):
    val = pulp.value(var)
    return val is not None and val > 0.5


def compatible_bases(symbol):
    return COMPAT.get(symbol, {symbol}) & set(BASE_ORDER)


def candidate_positions(n, start_seq, probe):
    """
    Return positions where the probe does not conflict with the known prefix.
    """
    positions = []
    k = len(probe)

    for start in range(n - k + 1):
        feasible = True

        for offset, symbol in enumerate(probe):
            pos = start + offset

            if pos < len(start_seq) and not match(start_seq[pos], symbol):
                feasible = False
                break

        if feasible:
            positions.append(start)

    return positions


def solve_exact_ilp(n, k, start_seq, probes, solver_name=None, threads=1):
    """
    Solve the SBH problem using an exact ILP model.

    b[p,a] = 1 if nucleotide a is selected at sequence position p.
    z[i,t] = 1 if probe i starts at sequence position t.
    y[i]   = 1 if probe i is covered by the resulting sequence.

    The model maximizes the number of distinct covered probes. Each selected
    probe placement is linked directly to the sequence nucleotides while
    fully respecting IUPAC compatibility. This avoids the graph-model error
    where pairwise probe compatibility did not guarantee one consistent
    concretization for all overlapping probes.
    """
    if not probes:
        return fallback_sequence(n, start_seq)

    if len(start_seq) >= n:
        return concretize_sequence(start_seq)[:n]

    model = pulp.LpProblem("Exact_Binary_Chip_Err_Minus_SBH", pulp.LpMaximize)
    bases = {
        (pos, base): pulp.LpVariable(f"b_{pos}_{base}", cat="Binary")
        for pos in range(n)
        for base in BASE_ORDER
    }

    for pos in range(n):
        model += pulp.lpSum(bases[(pos, base)] for base in BASE_ORDER) == 1

    concrete_start = concretize_sequence(start_seq)[:n]
    for pos, base in enumerate(concrete_start):
        model += bases[(pos, base)] == 1

    placements = {}
    covered = {}

    for i, probe in enumerate(probes):
        covered[i] = pulp.LpVariable(f"y_{i}", cat="Binary")
        positions = candidate_positions(n, concrete_start, probe)

        if not positions:
            model += covered[i] == 0
            continue

        probe_placements = []

        for start in positions:
            var = pulp.LpVariable(f"z_{i}_{start}", cat="Binary")
            placements[(i, start)] = var
            probe_placements.append(var)

            for offset, symbol in enumerate(probe):
                allowed = compatible_bases(symbol)
                model += var <= pulp.lpSum(
                    bases[(start + offset, base)] for base in allowed
                )

        model += pulp.lpSum(probe_placements) == covered[i]

    model += pulp.lpSum(covered.values())

    if solver_name == 'gurobi' and hasattr(pulp, 'GUROBI_CMD'):
        solver = pulp.GUROBI_CMD(msg=False)
    elif solver_name == 'cplex' and hasattr(pulp, 'CPLEX_CMD'):
        solver = pulp.CPLEX_CMD(msg=False)
    else:
        if TIME_LIMIT_SECONDS is None:
            try:
                solver = pulp.PULP_CBC_CMD(msg=False, threads=threads)
            except TypeError:
                solver = pulp.PULP_CBC_CMD(msg=False)
        else:
            try:
                solver = pulp.PULP_CBC_CMD(msg=False, timeLimit=TIME_LIMIT_SECONDS, threads=threads)
            except TypeError:
                solver = pulp.PULP_CBC_CMD(msg=False, timeLimit=TIME_LIMIT_SECONDS)

    status = model.solve(solver)
    try:
        status_name = pulp.LpStatus[status]
    except Exception:
        status_name = str(status)

    if DEBUG:
        sys.stderr.write(
            f"DEBUG: status={status_name} n={n} probes={len(probes)} "
            f"placements={len(placements)}\n"
        )

    has_solution = any(pulp.value(var) is not None for var in bases.values())
    if status_name not in {"Optimal", "Not Solved"} and not has_solution:
        return fallback_sequence(n, start_seq)

    dna = []
    for pos in range(n):
        chosen = next(
            (base for base in BASE_ORDER if selected_value(bases[(pos, base)])),
            'A',
        )
        dna.append(chosen)

    return ''.join(dna)


def main():
    data = sys.stdin.read()

    if not data.strip():
        return

    n, k, start_seq, probes = parse_xml(data)

    if DEBUG:
        sys.stderr.write(f"DEBUG main: n={n} k={k} start_seq={start_seq!r} probes_count={len(probes)}\n")

    result = solve_exact_ilp(n, k, start_seq, probes)

    sys.stdout.write(result + "\n")


if __name__ == "__main__":
    main()
