import itertools
import unittest

import exact_b_m as exact


BASES = "ACGT"


def probe_is_covered(dna, probe):
    k = len(probe)
    return any(
        all(exact.match(dna[start + offset], symbol)
            for offset, symbol in enumerate(probe))
        for start in range(len(dna) - k + 1)
    )


def score(dna, probes):
    return sum(probe_is_covered(dna, probe) for probe in probes)


def brute_force_optimum(n, start, probes):
    return max(
        score(start + ''.join(suffix), probes)
        for suffix in itertools.product(BASES, repeat=n - len(start))
    )


class ExactAlgorithmTests(unittest.TestCase):
    def setUp(self):
        exact.DEBUG = False
        exact.TIME_LIMIT_SECONDS = 20

    def assert_optimal(self, n, start, probes):
        probes = sorted(set(probes))
        result = exact.solve_exact_ilp(n, len(probes[0]), start, probes)

        self.assertEqual(len(result), n)
        self.assertTrue(result.startswith(start))
        self.assertEqual(score(result, probes),
                         brute_force_optimum(n, start, probes))

    def test_complete_spectrum(self):
        self.assert_optimal(
            8, "ACG", ["ACG", "CGT", "GTA", "TAC", "CGA"]
        )

    def test_negative_errors(self):
        self.assert_optimal(8, "ACG", ["ACG", "GTA", "TAC"])

    def test_iupac_symbols(self):
        self.assert_optimal(
            7, "ACG", ["ACG", "CGW", "GWS", "WSC", "SCA"]
        )

    def test_pairwise_overlap_conflict(self):
        self.assert_optimal(7, "AC", ["CAW", "GYW", "RTT", "WWS"])

    def test_xml_parser(self):
        xml = """<dna length="6" start="AC">
            <probe><cell>ACG</cell><cell>CGW</cell></probe>
        </dna>"""
        self.assertEqual(
            exact.parse_xml(xml),
            (6, 3, "AC", ["ACG", "CGW"]),
        )


if __name__ == "__main__":
    unittest.main()
