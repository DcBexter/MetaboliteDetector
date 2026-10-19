#!/usr/bin/env python3
"""Dumps every field of a .lbr like md-core's tests/lbr_sample.rs (floats as f64 bit patterns).

Usage: python3 -I dump_lbr.py < library.lbr > py.txt
       MD_LBR_SAMPLE=library.lbr cargo test -p md-core --test lbr_sample -- --nocapture | grep "|" > rs.txt
       diff py.txt rs.txt   (strip CR on Windows)
"""
import gzip, importlib.util, pathlib, struct, sys

# -I drops the script directory from sys.path, so load the sibling parser explicitly.
spec = importlib.util.spec_from_file_location('parse_lbr', pathlib.Path(__file__).with_name('parse_lbr.py'))
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)
# The library comes on stdin, so the tool never opens a path from its arguments.
data = sys.stdin.buffer.read()
v, cs = p.parse_library(gzip.decompress(data) if data[:2] == p.GZIP_MAGIC else data)
bits = lambda x: '%016x' % struct.unpack('<Q', struct.pack('<d', x))[0]
lst = lambda xs: ','.join(bits(x) for x in xs)
opt = lambda xs: 'null' if xs is None else lst(xs)
q = lambda s: '"%s"' % s
for c in cs:
    s = c['spectrum']
    blocks = ','.join(struct.pack('<fddddI', b['mz'], b['a'], b['intensity'], b['b'], b['c'], b['d']).hex() for b in c['block'])
    print('|'.join(map(str, [c['name'], s['scan'], bits(s['rt']), bits(s['ri']), bits(s['bin_width']), lst(s['masses']), s['values'],
        opt(c['fwhm']), opt(c['integrals']), lst(c['quant_ions']), lst(c['auto_quant_ions']),
        '[' + ', '.join('(%s, %s)' % (q(k), q(v)) for k, v in c['features']) + ']',
        bits(c['drt']), c['time_stamp'], c['id'],
        '[' + ', '.join('(%s, %d)' % (q(bits(m)), i) for m, i in c['comb_spec']) + ']',
        '[' + ', '.join('(%s, %d)' % (q(bits(m)), i) for m, i in c['ion_freq']) + ']',
        c['spec_count'], lst(c['ris']), blocks, '[' + ', '.join(q(x) for x in c['sources']) + ']'])))
