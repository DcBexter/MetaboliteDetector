# Cases for make-anova.sh: fixed edge cases and seeded random groups ("case NAME", "group V V ...", "end").
import random
random.seed(203)
cases = [
  ("hand", [[1, 2, 3], [2, 3, 4], [5, 6, 7]]),
  ("two_groups", [[1, 2, 3, 4], [2, 4, 6, 9]]),
  ("one_group", [[1, 2, 3]]),
  ("none", []),
  ("unequal", [[1, 2, 3], [2, 4], [5, 6, 7, 8], [3, 3, 4]]),
  ("zeros", [[0, 0, 5], [0, 0, 0], [0, 1, 0]]),
  ("separated", [[1, 1], [2, 2]]),
  ("equal", [[3, 3], [3, 3]]),
  ("equal_tenth", [[0.1, 0.1, 0.1], [0.1, 0.1]]),
  ("n_equals_k", [[1], [2], [3]]),
  ("empty_group", [[1, 2], []]),
  ("singletons_and_pair", [[1], [2, 3]]),
  ("large", [[1e6, 1e6 + 1, 1e6 + 2], [1e6 + 1, 1e6 + 2, 1e6 + 3]]),
  ("tiny_p", [[1, 1.01, 0.99, 1.02], [5, 5.01, 4.99, 5.02], [9, 9.01, 8.99, 9.02]]),
]
for i in range(60):
    k = random.randint(2, 6)
    shift = random.choice([0, 0, 1, 5])
    groups = [[round(random.gauss(100 + shift * g, 20), 3) for _ in range(random.randint(2, 8))] for g in range(k)]
    cases.append((f"random{i}", groups))
for name, groups in cases:
    print("case", name)
    for g in groups:
        print("group", *g)
    print("end")
