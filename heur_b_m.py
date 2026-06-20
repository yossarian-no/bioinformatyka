#!/usr/bin/env python3

import sys
import time
import random
import xml.etree.ElementTree as ET


COMPAT = {
    'A': {'A'}, 'T': {'T'}, 'G': {'G'}, 'C': {'C'},

    'R': {'A', 'G'},
    'Y': {'C', 'T'},

    'S': {'G', 'C'},
    'W': {'A', 'T'},
}

BASE_ORDER = ['A', 'C', 'G', 'T']
RNG_SEED = 156202
TIME_LIMIT_SECONDS = 95.0


def match(a, b):
    """
    Checks whether two symbols are compatible.
    Works for ordinary DNA symbols and binary/IUPAC symbols.
    """
    return bool(COMPAT.get(a, {a}) & COMPAT.get(b, {b}))


def concretize_symbol(ch):
    """
    Converts ambiguous symbol to one concrete nucleotide.
    Deterministic choice according to BASE_ORDER.
    """
    options = COMPAT.get(ch, {ch})

    for b in BASE_ORDER:
        if b in options:
            return b

    return 'A'


def concretize_sequence(seq):
    """
    Converts a sequence possibly containing ambiguous symbols
    into a concrete A/C/G/T sequence.
    """
    return ''.join(concretize_symbol(ch) for ch in seq)


def parse_xml(data):
    """
    Reads XML from stdin.

    Expected format:
        <dna length="500" start="...">
            <probe pattern="...">
                <cell>...</cell>
                ...
            </probe>
        </dna>

    Only <cell> values are treated as probes.
    The attribute pattern is ignored.
    """
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
        n = len(start_seq) if start_seq else k

    return n, k, start_seq, probes


def min_shift_from_sequence(dna, probe):
    """
    Finds minimal shift needed to append probe after current DNA sequence.

    cost = shift = number of new nucleotides added to DNA

    For example:
        dna suffix: ACGT
        probe:        GTA
        overlap: GT
        cost: 1
    """
    k = len(probe)
    max_overlap = min(len(dna), k - 1)

    for shift in range(1, k + 1):
        overlap_len = k - shift

        if overlap_len > max_overlap:
            continue

        if overlap_len <= 0:
            return k

        start = len(dna) - overlap_len

        ok = True
        for i in range(overlap_len):
            if not match(dna[start + i], probe[i]):
                ok = False
                break

        if ok:
            return shift

    return k


def min_shift_between_probes(a, b):
    """
    Minimal shift for placing probe b after probe a.
    Used only for heuristic scoring.
    """
    k = len(b)

    for shift in range(1, k + 1):
        overlap_len = k - shift

        if overlap_len <= 0:
            return k

        if overlap_len > len(a):
            continue

        ok = True
        start = len(a) - overlap_len

        for i in range(overlap_len):
            if not match(a[start + i], b[i]):
                ok = False
                break

        if ok:
            return shift

    return k


def append_probe(dna, probe, cost, n):
    """
    Appends only the new part of a probe to current DNA.
    """
    remaining = n - len(dna)

    if remaining <= 0:
        return dna

    cost = min(cost, remaining)

    if cost <= 0:
        return dna

    addition = ''.join(concretize_symbol(ch) for ch in probe[-cost:])

    return dna + addition


def probe_covered_at(dna, pos, probe):
    """
    Checks whether probe is compatible with DNA substring at position pos.
    """
    k = len(probe)

    if pos < 0 or pos + k > len(dna):
        return False

    for i in range(k):
        if not match(dna[pos + i], probe[i]):
            return False

    return True


def probe_covered_anywhere(dna, probe):
    """
    Checks whether probe occurs in DNA in compatible form.
    """
    k = len(probe)

    if len(dna) < k:
        return False

    for pos in range(0, len(dna) - k + 1):
        if probe_covered_at(dna, pos, probe):
            return True

    return False


def update_covered_after_append(dna, probes, covered):
    """
    Updates the set of covered probes.

    Since DNA is built by appending characters at the end, most new probes
    can only appear near the end. This is faster than scanning whole DNA
    after every step.
    """
    if not probes:
        return covered

    k = len(probes[0])
    start_pos = max(0, len(dna) - k)

    for idx, probe in enumerate(probes):
        if idx in covered:
            continue

        found = False

        for pos in range(start_pos, len(dna)):
            if probe_covered_at(dna, pos, probe):
                found = True
                break

        if found:
            covered.add(idx)

    return covered


def initial_covered(dna, probes):
    """
    Finds probes already covered by the initial DNA fragment.
    """
    covered = set()

    for idx, probe in enumerate(probes):
        if probe_covered_anywhere(dna, probe):
            covered.add(idx)

    return covered


def count_covered(dna, probes):
    """
    Counts how many probes are covered by the final DNA sequence.
    Used to compare several heuristic runs.
    """
    score = 0

    for probe in probes:
        if probe_covered_anywhere(dna, probe):
            score += 1

    return score


def compute_future_score(probes):
    """
    Precomputes a simple heuristic score for each probe:
    how many other probes can follow it with a small shift.

    This does not solve the problem exactly. It is only used to break ties
    in the greedy algorithm.
    """
    m = len(probes)

    if m == 0:
        return []

    k = len(probes[0])
    good_limit = min(4, k)

    future = [0] * m

    for i in range(m):
        cnt = 0

        for j in range(m):
            if i == j:
                continue

            c = min_shift_between_probes(probes[i], probes[j])

            if c <= good_limit:
                cnt += 1

        future[i] = cnt

    return future


def choose_start_probe(probes, future_score):
    """
    Chooses a start probe if XML does not contain start sequence.
    """
    if not probes:
        return None

    best_idx = 0
    best_score = future_score[0] if future_score else 0

    for i in range(1, len(probes)):
        sc = future_score[i] if future_score else 0

        if sc > best_score:
            best_score = sc
            best_idx = i

    return best_idx


def greedy_construct(probes, n, start_seq, future_score, rng, randomized=False):
    """
    Greedy construction of DNA.

    At every step, choose a not-yet-covered probe that can be appended
    with the smallest cost. Ties are broken by future_score and optionally
    by randomization.
    """
    if n <= 0:
        return ''

    if start_seq:
        dna = concretize_sequence(start_seq)[:n]
        covered = initial_covered(dna, probes)
    else:
        start_idx = choose_start_probe(probes, future_score)

        if start_idx is None:
            dna = 'A' * n
            covered = set()
        else:
            dna = concretize_sequence(probes[start_idx])[:n]
            covered = initial_covered(dna, probes)

    if len(dna) >= n:
        return dna[:n]

    m = len(probes)

    while len(dna) < n:
        remaining = n - len(dna)
        candidates = []

        for idx, probe in enumerate(probes):
            if idx in covered:
                continue

            cost = min_shift_from_sequence(dna, probe)

            if cost <= remaining:
                fs = future_score[idx] if future_score else 0
                candidates.append((cost, -fs, idx))

        if not candidates:
            dna += 'A' * remaining
            break

        candidates.sort()

        if randomized and len(candidates) > 1:
            pool_size = min(8, len(candidates))
            chosen = rng.choice(candidates[:pool_size])
        else:
            chosen = candidates[0]

        cost, _, idx = chosen

        dna = append_probe(dna, probes[idx], cost, n)
        covered.add(idx)
        covered = update_covered_after_append(dna, probes, covered)

    return dna[:n]


def improve_by_point_mutations(
    dna, probes, n, rng, deadline, fixed_prefix_length=0
):
    """
    Small local improvement.

    It tries single nucleotide substitutions and keeps a mutation only
    if it does not decrease the number of covered probes.

    This is heuristic and time-limited.
    """
    if len(dna) != n:
        dna = dna[:n].ljust(n, 'A')

    best = list(dna)
    best_score = count_covered(''.join(best), probes)

    if best_score == len(probes):
        return ''.join(best)

    fixed_prefix_length = min(max(fixed_prefix_length, 0), n)

    if fixed_prefix_length == n:
        return ''.join(best)

    attempts = min(300, max(30, n))

    for _ in range(attempts):
        if time.time() >= deadline:
            break

        pos = rng.randrange(fixed_prefix_length, n)
        old = best[pos]

        bases = BASE_ORDER[:]
        rng.shuffle(bases)

        for b in bases:
            if b == old:
                continue

            best[pos] = b
            candidate = ''.join(best)
            sc = count_covered(candidate, probes)

            if sc >= best_score:
                best_score = sc
                break

            best[pos] = old

    return ''.join(best)


def solve_heuristic(n, k, start_seq, probes):
    """
    Main heuristic solver:
    - deterministic greedy run,
    - several randomized greedy restarts,
    - optional local improvement.
    """
    if n <= 0:
        return ''

    if not probes:
        if start_seq:
            return concretize_sequence(start_seq)[:n].ljust(n, 'A')
        return 'A' * n

    rng = random.Random(RNG_SEED)
    start_time = time.time()
    deadline = start_time + TIME_LIMIT_SECONDS

    future_score = compute_future_score(probes)

    m = len(probes)

    if m > 1500:
        restarts = 3
    elif m > 800:
        restarts = 6
    elif m > 300:
        restarts = 10
    else:
        restarts = 20

    best_dna = None
    best_score = -1

    for r in range(restarts):
        if time.time() >= deadline:
            break

        randomized = r > 0

        dna = greedy_construct(
            probes=probes,
            n=n,
            start_seq=start_seq,
            future_score=future_score,
            rng=rng,
            randomized=randomized
        )

        if time.time() < deadline:
            dna = improve_by_point_mutations(
                dna,
                probes,
                n,
                rng,
                deadline,
                fixed_prefix_length=len(start_seq),
            )

        sc = count_covered(dna, probes)

        if sc > best_score:
            best_score = sc
            best_dna = dna

        if best_score == len(probes):
            break

    if best_dna is None:
        if start_seq:
            best_dna = concretize_sequence(start_seq)[:n].ljust(n, 'A')
        else:
            best_dna = 'A' * n

    if len(best_dna) < n:
        best_dna += 'A' * (n - len(best_dna))

    return best_dna[:n]


def main():
    data = sys.stdin.read()

    if not data.strip():
        return

    n, k, start_seq, probes = parse_xml(data)

    result = solve_heuristic(n, k, start_seq, probes)

    sys.stdout.write(result + '\n')


if __name__ == '__main__':
    main()
