# Cases for make-anova.sh: fixed edge cases and seeded random groups ("case NAME", "group V V ...", "end").

state = 203


def uniform():
    """A fixed linear congruential sequence in [0, 1): the cases must be the same on every run."""
    global state
    state = (state * 6364136223846793005 + 1442695040888963407) % 2**64
    return (state >> 11) / 2**53


def below(n):
    return int(uniform() * n)


def normal(mean, sd):
    return mean + sd * (sum(uniform() for _ in range(12)) - 6)


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
    k = 2 + below(5)
    shift = [0, 0, 1, 5][below(4)]
    groups = [[round(normal(100 + shift * g, 20), 3) for _ in range(2 + below(7))] for g in range(k)]
    cases.append((f"random{i}", groups))
for name, groups in cases:
    print("case", name)
    for g in groups:
        print("group", *g)
    print("end")
