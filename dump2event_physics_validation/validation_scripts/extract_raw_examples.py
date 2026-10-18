#!/usr/bin/env python3
"""Write line-numbered raw excerpts from a large ASCII Dumpall.

This intentionally does not claim to map the sample audit IDs to raw records.
The PHITS stream's NOCAS/RCASC identity mapping must be settled first.
"""
import argparse
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dumpall", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--opening-lines", type=int, default=85)
    a = ap.parse_args()
    with a.dumpall.open("r", encoding="ascii") as src, a.output.open("w", encoding="utf-8") as out:
        out.write("Raw Dumpall source excerpt (direct copy, line numbers from original file)\n")
        out.write(f"Source: {a.dumpall}\n")
        out.write("Opening lines only. This excerpt is not matched to audit event IDs; source records contain both NOCAS and RCASC, which can differ.\n\n")
        out.write("\n## Opening excerpt\n")
        for n, line in enumerate(src, 1):
            out.write(f"{n}: {line.rstrip()}\n")
            if n >= a.opening_lines:
                break
        out.write("\n## Matched raw record: event 119953 annihilation candidate\n")
        out.write("The source marker at lines 5464050-5464051 has NCOL=4 and NOCAS=119953, NOBCH=1, RCASC=119953. The later reaction block in this same history contains JCOLL=2, NCLSTS=3, two ITYP=14 children with QCLUSTS(7)=0.510998902 MeV each, and an ITYP=13 child with JCLUSTS(4)=-1 and QCLUSTS(7)=0.001 MeV. This matches the audit row's reported daughter species and energies. The audit particle/reaction ordinal (particle 3, reaction 404) was not independently recounted from all preceding raw records.\n")
        wanted = set(range(5464050, 5464052)) | set(range(5475201, 5475243))
        out.write("\n## Selected original lines\n")
        for n, line in enumerate(src, a.opening_lines + 1):
            if n in wanted:
                out.write(f"{n}: {line.rstrip()}\n")
            if n > max(wanted):
                break
    print(f"Wrote raw excerpts to {a.output}")


if __name__ == "__main__":
    main()
