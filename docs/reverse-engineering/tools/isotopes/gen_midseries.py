# Synthetic cases for make-midseries.sh: random integrals for fragments of
# library substance 0. Runs listed in GAPS have no integrals for that case (NaN in the old program).
import random
random.seed(28)
cases = [
  # name, compound, tracer, flag, runs, frags [(lo, hi, formula, size)]
  ("lib_c", 0, "C", 1, 3, [(57, 60, "-", 0), (83, 86, "-", 0)]),
  ("lib_n", 0, "N", 1, 2, [(71, 74, "-", 0)]),
  ("lib_rev", 0, "C", 1, 2, [(99, 97, "-", 0)]),
  ("form_c", 0, "C", 0, 3, [(57, 61, "C4H9", 0), (85, 88, "C5H9O", 3)]),
  ("form_none", 0, "-", 0, 2, [(71, 74, "C4H7O", 2)]),
  ("form_big", 0, "C", 0, 2, [(57, 59, "C4H9", 6)]),
  ("form_gap", 0, "C", 0, 4, [(71, 74, "C4H7O", 3)]),
  ("lib_one", 0, "C", 1, 2, [(57, 59, "-", 0)]),
]
GAPS = {"form_gap": {1}, "lib_one": {0}}
out = []
for name, comp, tracer, flag, runs, frags in cases:
    out.append(f"case {name} {comp} {tracer} {flag} {runs}")
    for lo, hi, f, s in frags:
        out.append(f"frag {lo} {hi} {f} {s}")
    for lo, hi, f, s in frags:
        a, b = min(lo, hi), max(lo, hi)
        for r in range(runs):
            base = random.uniform(1e4, 1e6)
            for k, ion in enumerate(range(a, b + 1)):
                v = round(base * random.uniform(0.05, 1.0) / (k + 1), 3)
                if r not in GAPS.get(name, ()):
                    out.append(f"int {r} {ion} {v}")
    out.append("end")
print("\n".join(out))
