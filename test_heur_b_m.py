import unittest

import heur_b_m as heuristic


class HeuristicAlgorithmTests(unittest.TestCase):
    def test_known_start_is_preserved(self):
        start = "ACG"
        probes = ["ACG", "CGT", "GTA", "TAC", "CGA"]

        result = heuristic.solve_heuristic(8, 3, start, probes)

        self.assertEqual(len(result), 8)
        self.assertTrue(result.startswith(start))


if __name__ == "__main__":
    unittest.main()
