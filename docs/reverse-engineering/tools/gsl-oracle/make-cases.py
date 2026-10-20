#!/usr/bin/env python3
"""Black-box cases of the old quadratic fit for the ignored test quadfit::matches_the_old_fit_bit_for_bit.

  docs/reverse-engineering/tools/gsl-oracle/make-cases.py [N] > cases.txt
  MD_QUADFIT_CASES=cases.txt cargo test -p md-core quadfit -- --ignored

Feeds N random precise-apex fits (scan, bins, three intensities) through run.sh and prints one case per
line: x0 x1 x2 y0 y1 y2 c0 c1 c2 as IEEE bit patterns in hex. The outputs come from the old program:
keep them local, do not commit them.
"""
import os
import random
import struct
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
bits = lambda v: '%016x' % struct.unpack('<Q', struct.pack('<d', v))[0]
f32 = lambda v: struct.unpack('<f', struct.pack('<f', v))[0]
rng = random.Random(250)  # reproducible test cases, not security
n = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
inputs = []
for _ in range(n):
    bins = rng.choice([10, 10, 10, rng.randint(1, 100)])
    scan = rng.choice([rng.randint(1, 100), rng.randint(1, 10000), rng.randint(1, 1000000)])
    k = rng.random()
    if k < 0.5:  # integer intensities, often near-ties
        y = [float(rng.randint(0, 60000)) for _ in range(3)]
    elif k < 0.8:  # f32 intensities, also negative like matched-filter output
        y = [f32(rng.uniform(-1e4, 1e5)) for _ in range(3)]
    else:  # flat tops
        y = [float(rng.randint(0, 1000))] * 3
    inputs.append([float(bins * (scan - 1 + i)) for i in range(3)] + y)
text = ''.join('fit ' + ' '.join(v.hex() for v in a) + '\n' for a in inputs)
out = subprocess.run([os.path.join(here, 'run.sh')], input=text, capture_output=True, text=True, check=True).stdout
lines = out.splitlines()
assert len(lines) == n, (len(lines), n)
for a, line in zip(inputs, lines):
    print(' '.join([bits(v) for v in a] + line.split()[1:]))
