#!/usr/bin/env python3
"""Parser for MetaboliteDetector 3.3 compound libraries (.lbr).

The file is a gzip stream. The uncompressed payload is little-endian:

    uint32 version            (sample file: 2; shipped writer emits 5, see file-formats.md)
    uint32 compound_count
    compound_count x LibraryCompound records

Only the version-2 record layout observed in the sample library is implemented.
Version 5 adds a Peak<double,short> model block and trailing formulas/mid_sizes
maps, which are not decoded here.

Strings are uint32 length + raw bytes (no terminator). They are decoded as
latin-1 so the byte count is preserved; the real encoding is unverified.

Usage (run from any directory, no local imports needed):
    python3 -I parse_lbr.py <library.lbr | decompressed.raw>
"""

import gzip
import struct
import sys

GZIP_MAGIC = b"\x1f\x8b"
SUPPORTED_VERSIONS = (2,)
MASK32 = 0xFFFFFFFF


class Reader:
    def __init__(self, data):
        self.data = data
        self.pos = 0

    def take(self, n):
        end = self.pos + n
        if n < 0 or end > len(self.data):
            raise ValueError(f"truncated record: need {n} bytes at offset {self.pos}")
        chunk = self.data[self.pos:end]
        self.pos = end
        return chunk

    def unpack(self, fmt):
        return struct.unpack(fmt, self.take(struct.calcsize(fmt)))[0]

    def u32(self):
        return self.unpack("<I")

    def i32(self):
        return self.unpack("<i")

    def i64(self):
        return self.unpack("<q")

    def f32(self):
        return self.unpack("<f")

    def f64(self):
        return self.unpack("<d")

    def string(self):
        return self.take(self.u32()).decode("latin-1")


def read_f64_vector_or_null(r):
    """int32 count (-1 = null pointer) followed by count float64 values."""
    n = r.i32()
    if n < 0:
        return None
    return [r.f64() for _ in range(n)]


def read_f32_set(r):
    """uint32 count followed by count float32 values (std::set<float>, sorted)."""
    return [r.f32() for _ in range(r.u32())]


def read_spectrum(r):
    scan = r.u32()
    rt = r.f64()
    ri = r.f64()
    bin_width = r.f32()
    n_masses = r.u32()
    masses = [r.f32() for _ in range(n_masses)]
    n_values = r.u32()
    values = [r.i32() for _ in range(n_values)]
    check = r.u32()
    if check != (scan + n_masses) & MASK32:
        raise ValueError(f"spectrum check mismatch: {check} != scan+n_masses")
    return {
        "scan": scan,
        "rt": rt,
        "ri": ri,
        "bin_width": bin_width,
        "masses": masses,
        "values": values,
    }


def read_compound(r):
    n_peaks = r.u32()
    if n_peaks:
        raise NotImplementedError("PeakEssential<float> records are not decoded")
    name = r.string()
    fwhm = read_f64_vector_or_null(r)
    integrals = read_f64_vector_or_null(r)
    quant_ions = read_f32_set(r)
    auto_quant_ions = read_f32_set(r)
    features = [(r.string(), r.string()) for _ in range(r.u32())]
    drt = r.f64()
    time_stamp = r.i64()
    ident = r.i64()
    check = r.u32()
    expected = len(quant_ions) + len(features) + len(name.encode("latin-1")) + n_peaks
    if check != expected & MASK32:
        raise ValueError(f"compound check mismatch for {name!r}: {check} != {expected}")
    return {
        "name": name,
        "spectrum": None,  # filled by caller
        "fwhm": fwhm,
        "integrals": integrals,
        "quant_ions": quant_ions,
        "auto_quant_ions": auto_quant_ions,
        "features": features,
        "drt": drt,
        "time_stamp": time_stamp,
        "id": ident,
    }


def read_block_record(r):
    """40-byte record of unknown meaning: float mz, f64 1.0, f64 intensity, f64 0.0, f64 -1.0, u32 1.

    The intensity equals the comb_spec intensity for the same m/z (see file-formats.md, Unknowns).
    """
    mz = r.f32()
    a = r.f64()
    intensity = r.f64()
    b = r.f64()
    c = r.f64()
    d = r.u32()
    return {"mz": mz, "a": a, "intensity": intensity, "b": b, "c": c, "d": d}


def read_library_compound(r):
    spectrum = read_spectrum(r)
    compound = read_compound(r)
    compound["spectrum"] = spectrum
    compound["comb_spec"] = [(r.f32(), r.i32()) for _ in range(r.u32())]
    compound["ion_freq"] = [(r.f32(), r.u32()) for _ in range(r.u32())]
    compound["spec_count"] = r.u32()
    compound["ris"] = [r.f64() for _ in range(r.u32())]
    compound["block"] = [read_block_record(r) for _ in range(r.u32())]
    compound["sources"] = [r.string() for _ in range(r.u32())]
    return compound


def parse_library(data):
    r = Reader(data)
    version = r.u32()
    if version not in SUPPORTED_VERSIONS:
        raise NotImplementedError(f"library version {version} layout not decoded")
    count = r.u32()
    compounds = [read_library_compound(r) for _ in range(count)]
    if r.pos != len(data):
        raise ValueError(f"{len(data) - r.pos} trailing bytes after last compound")
    return version, compounds


def describe(c):
    sp = c["spectrum"]
    return (
        f"name={c['name']!r} ri={sp['ri']:.3f} rt_field={sp['rt']:.6g} "
        f"scan={sp['scan']} masses={len(sp['masses'])} comb_spec={len(c['comb_spec'])} "
        f"ion_freq={len(c['ion_freq'])} ris={c['ris']} quant_ions={c['quant_ions']} "
        f"features={len(c['features'])} sources={c['sources']}"
    )


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    with open(argv[1], "rb") as fh:
        data = fh.read()
    if data[:2] == GZIP_MAGIC:
        data = gzip.decompress(data)
    version, compounds = parse_library(data)
    print(f"version={version} compounds={len(compounds)} bytes={len(data)} (fully consumed)")
    for c in compounds[:3]:
        print(describe(c))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
