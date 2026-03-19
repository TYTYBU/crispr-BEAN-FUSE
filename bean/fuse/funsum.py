"""
FUNSUM matrix loading utilities.

FUNSUM (FUNctional SUbstitution Matrix) encodes the expected effect of each
amino-acid substitution.  Two flavours are bundled with the package:

* ``funsum``        – 21×21 matrix including stop-gain (*) variants
* ``funsum_noLOF``  – 20×20 matrix excluding stop-gain variants
* ``funsum_SS``     – per-secondary-structure dict of the above (with *)
* ``funsum_SS_noLOF`` – per-secondary-structure dict (without *)
"""

from __future__ import annotations

import os
import pickle
from typing import Dict

import pandas as pd

_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

_CACHE: Dict[str, object] = {}


def _load_pkl(name: str):
    if name not in _CACHE:
        path = os.path.join(_DATA_DIR, f"{name}.pkl")
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"Bundled FUNSUM file not found: {path}. "
                "Re-run the conversion step to regenerate package data."
            )
        with open(path, "rb") as fh:
            _CACHE[name] = pickle.load(fh)
    return _CACHE[name]


def load_funsum(include_lof: bool = True) -> pd.DataFrame:
    """Return the overall FUNSUM matrix as a DataFrame.

    Parameters
    ----------
    include_lof:
        If *True* (default) return the 21×21 matrix that includes stop-gain
        variants.  If *False* return the 20×20 matrix without them.
    """
    key = "funsum" if include_lof else "funsum_noLOF"
    return _load_pkl(key).copy()


def load_funsum_ss(include_lof: bool = True) -> Dict[str, pd.DataFrame]:
    """Return the per-secondary-structure FUNSUM dictionary.

    Keys are DSSP single-character codes (G, H, I, E, B, T, S, C) plus the
    grouped labels Helices, Strands, Loops.

    Parameters
    ----------
    include_lof:
        Same semantics as :func:`load_funsum`.
    """
    key = "funsum_SS" if include_lof else "funsum_SS_noLOF"
    raw = _load_pkl(key)
    return {k: v.copy() for k, v in raw.items()}


def funsum_to_sub_table(df_funsum: pd.DataFrame) -> pd.DataFrame:
    """Flatten a FUNSUM matrix into a long-format substitution table.

    Returns a DataFrame with columns ``aa_pair`` and ``score``, where
    ``aa_pair`` is the two-character string ``ref + alt``.
    """
    rows = []
    for ref in df_funsum.index:
        for alt in df_funsum.columns:
            rows.append({"aa_pair": f"{ref}{alt}", "score": df_funsum.loc[ref, alt]})
    return pd.DataFrame(rows)
