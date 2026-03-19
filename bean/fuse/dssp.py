"""
DSSP .dss file parser.

Parses the fixed-width residue table produced by the DSSP program
(Kabsch & Sander, Biopolymers 1983) into a tidy DataFrame.

Only the columns needed for FUSE scoring are extracted:
  - ``resnum``  : sequential residue number (1-based, matches aapos in BEAN)
  - ``respdb``  : PDB residue sequence number (may differ from resnum)
  - ``chain``   : chain identifier
  - ``aa``      : one-letter amino acid code
  - ``ss``      : DSSP secondary structure code (H, B, E, G, I, T, S, or blank)
  - ``sasa``    : solvent accessible surface area (Å²)
"""

from __future__ import annotations

import re
from typing import Optional

import pandas as pd


# The residue table starts after the header line that begins with "  #  RESIDUE"
_HEADER_PATTERN = re.compile(r"^\s+#\s+RESIDUE")

# DSSP fixed-width column positions (0-based, exclusive end)
# Reference: DSSP format specification
_COL_RESNUM = (0, 5)       # sequential residue number
_COL_RESPDB = (5, 10)      # PDB residue number
_COL_CHAIN  = (11, 12)     # chain
_COL_AA     = (13, 14)     # one-letter amino acid (! = chain break)
_COL_SS     = (16, 17)     # secondary structure
_COL_SASA   = (35, 38)     # solvent accessibility (integer)


def parse_dssp(path: str) -> pd.DataFrame:
    """Parse a DSSP ``.dss`` file and return a tidy residue DataFrame.

    Parameters
    ----------
    path:
        Path to the DSSP output file (extension ``.dss`` or ``.dssp``).

    Returns
    -------
    pd.DataFrame
        Columns: ``resnum``, ``respdb``, ``chain``, ``aa``, ``ss``, ``sasa``.
        ``ss`` is the single DSSP code; blank (' ') is stored as ``'C'``
        (coil/loop) to be consistent with the convention used in FUSE.
        Chain-break records (``aa == '!'``) are dropped.
    """
    with open(path, "r") as fh:
        lines = fh.readlines()

    # Find the header line
    header_idx: Optional[int] = None
    for i, line in enumerate(lines):
        if _HEADER_PATTERN.match(line):
            header_idx = i
            break

    if header_idx is None:
        raise ValueError(
            f"Could not locate the residue table header in DSSP file: {path}"
        )

    records = []
    for line in lines[header_idx + 1 :]:
        if len(line) < _COL_SASA[1]:
            continue  # skip short / blank lines

        aa = line[_COL_AA[0] : _COL_AA[1]].strip()
        if not aa or aa == "!":
            # chain-break record – skip
            continue

        resnum_s = line[_COL_RESNUM[0] : _COL_RESNUM[1]].strip()
        respdb_s = line[_COL_RESPDB[0] : _COL_RESPDB[1]].strip()
        chain    = line[_COL_CHAIN[0]  : _COL_CHAIN[1] ].strip()
        ss       = line[_COL_SS[0]     : _COL_SS[1]    ].strip()
        sasa_s   = line[_COL_SASA[0]  : _COL_SASA[1]  ].strip()

        if not ss:
            ss = "C"  # blank = irregular / coil

        try:
            resnum = int(resnum_s)
            respdb = int(re.sub(r"[A-Za-z]", "", respdb_s)) if respdb_s else resnum
            sasa   = float(sasa_s) if sasa_s else float("nan")
        except ValueError:
            continue

        records.append(
            {
                "resnum": resnum,
                "respdb": respdb,
                "chain":  chain,
                "aa":     aa,
                "ss":     ss,
                "sasa":   sasa,
            }
        )

    if not records:
        raise ValueError(f"No residue records parsed from DSSP file: {path}")

    return pd.DataFrame(records)
