"""
Utility functions for FUSE score computation.

Includes:
  - Amino-acid code conversion (three-letter ↔ one-letter)
  - Edit-string parsing from bean element result CSVs
  - Functional-class annotation (SYN / MIS / LOF)
  - Score normalisation (mirrors normalize_scoreset() in the original R code)
  - James-Stein shrinkage estimator (mirrors get_js() in the original R code)
"""

from __future__ import annotations

import re
import warnings
from typing import Optional

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Amino-acid code tables
# ---------------------------------------------------------------------------

THREE_TO_ONE: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
    # stop codon aliases sometimes seen in DMS data
    "TER": "*", "STOP": "*", "TRM": "*",
}

ONE_LETTER_SET = set("ACDEFGHIKLMNPQRSTVWY*")


def to_one_letter(code: str) -> str:
    """Convert an amino-acid code to its single-letter representation.

    Accepts one-letter codes (returned unchanged), three-letter codes
    (case-insensitive), and ``'*'`` / ``'Z'`` for stop codons.

    Raises ``ValueError`` for unrecognised codes.
    """
    code = code.strip()
    if code == "Z":
        return "*"
    if len(code) == 1:
        upper = code.upper()
        if upper in ONE_LETTER_SET or upper == "*":
            return upper
        raise ValueError(f"Unrecognised single-letter amino acid code: '{code}'")
    if len(code) == 3:
        upper = code.upper()
        if upper in THREE_TO_ONE:
            return THREE_TO_ONE[upper]
        raise ValueError(f"Unrecognised three-letter amino acid code: '{code}'")
    raise ValueError(f"Cannot interpret amino acid code: '{code}'")


# ---------------------------------------------------------------------------
# Edit-string parsing
# ---------------------------------------------------------------------------

# Patterns for bean AminoAcidEdit string format:
#   "pos:ref>alt"          e.g. "123:A>G"
#   "gene:pos:ref>alt"     e.g. "LDLR:123:A>G"
# Also accept the pipe-separated CodingNoncodingAllele format:
#   "LDLR:123:A>G|..."    (aa part before the pipe is used)

_AA_EDIT_RE = re.compile(
    r"^(?:[^:]+:)?(-?\d+):([A-Za-z*]{1,3})>([A-Za-z*Z]{1,3})$"
)


def parse_edit_string(edit_str: str) -> tuple[int, str, str] | None:
    """Parse a bean amino-acid edit string into (aapos, aaref, aaalt).

    Handles formats:
      - ``"pos:ref>alt"``
      - ``"gene:pos:ref>alt"``
      - ``"gene:pos:ref>alt|nt_edit_string"``  (CodingNoncodingAllele)

    Returns ``None`` if the string cannot be parsed or ``aapos <= 0``.
    """
    if not isinstance(edit_str, str):
        return None

    # Strip the nucleotide part if present (CodingNoncodingAllele format)
    aa_part = edit_str.split("|")[0].strip()

    m = _AA_EDIT_RE.match(aa_part)
    if not m:
        return None

    pos_s, ref_s, alt_s = m.group(1), m.group(2), m.group(3)

    try:
        aapos = int(pos_s)
        aaref = to_one_letter(ref_s)
        aaalt = to_one_letter(alt_s)
    except (ValueError, KeyError):
        return None

    if aapos <= 0:
        return None

    return aapos, aaref, aaalt


def extract_aa_columns(
    df: pd.DataFrame,
    edit_col: str,
) -> pd.DataFrame:
    """Parse the edit column and add ``aapos``, ``aaref``, ``aaalt`` columns.

    Rows whose edit string cannot be parsed receive NaN / NaT values and are
    flagged with a warning.

    Parameters
    ----------
    df:
        Input DataFrame (typically a bean element result CSV).
    edit_col:
        Name of the column containing amino-acid edit strings.

    Returns
    -------
    pd.DataFrame
        Copy of *df* with three new columns added: ``aapos`` (int),
        ``aaref`` (str), ``aaalt`` (str).
    """
    df = df.copy()
    parsed = df[edit_col].map(parse_edit_string)

    n_failed = parsed.isna().sum()
    if n_failed > 0:
        warnings.warn(
            f"{n_failed} rows could not be parsed from column '{edit_col}' "
            "and will be excluded from FUSE scoring.",
            stacklevel=2,
        )

    df["aapos"] = parsed.map(lambda x: x[0] if x is not None else np.nan)
    df["aaref"] = parsed.map(lambda x: x[1] if x is not None else np.nan)
    df["aaalt"] = parsed.map(lambda x: x[2] if x is not None else np.nan)
    df["aapos"] = pd.to_numeric(df["aapos"], errors="coerce").astype("Int64")

    return df


# ---------------------------------------------------------------------------
# Functional class annotation
# ---------------------------------------------------------------------------

def annotate_functional_class(df: pd.DataFrame) -> pd.DataFrame:
    """Add a ``functional_class`` column (SYN / MIS / LOF).

    Requires columns ``aaref`` and ``aaalt``.
    """
    df = df.copy()
    df["functional_class"] = np.nan

    mis_mask = (df["aaref"] != df["aaalt"]) & (df["aaalt"] != "*")
    lof_mask = (df["aaref"] != "*") & (df["aaalt"] == "*")
    syn_mask = df["aaref"] == df["aaalt"]

    df.loc[mis_mask, "functional_class"] = "MIS"
    df.loc[lof_mask, "functional_class"] = "LOF"
    df.loc[syn_mask, "functional_class"] = "SYN"

    return df


# ---------------------------------------------------------------------------
# Score normalisation
# ---------------------------------------------------------------------------

def normalize_scoreset(
    df: pd.DataFrame,
    lower_bound: float = 0.1,
    upper_bound: float = 0.9,
    force_quantile: bool = True,
) -> pd.DataFrame:
    """Normalise raw scores to a [LOF=1, SYN=0]-anchored scale.

    Mirrors ``normalize_scoreset()`` from ``AA_proj_functions_101421.R``.

    The function:
      1. Determines the orientation of the score (higher = more damaging or
         less damaging) using the relative positions of the LOF and SYN medians
         (or falls back to Proline as a reference when both groups are absent).
      2. Subtracts the SYN median and divides by (LOF median − SYN median) so
         that SYN variants cluster around 0 and LOF variants cluster around 1.
      3. If ``force_quantile`` is *True*, re-normalises a second time using the
         lower/upper quantiles of the missense distribution (useful when the
         initial scale is already close to the target).

    Parameters
    ----------
    df:
        DataFrame with columns ``aaref``, ``aaalt``, ``functional_class``,
        and ``raw_score``.
    lower_bound, upper_bound:
        Quantiles of the missense distribution used as surrogate SYN / LOF
        anchors when one group is absent (or for ``force_quantile``).
    force_quantile:
        Re-normalise using missense quantiles after the primary normalisation.

    Returns
    -------
    pd.DataFrame
        Copy of *df* with a new ``norm_raw_score`` column.
    """
    df = df.copy()

    ind_lof = df["functional_class"] == "LOF"
    ind_mis = df["functional_class"] == "MIS"
    ind_syn = df["functional_class"] == "SYN"

    median_lof = df.loc[ind_lof, "raw_score"].median() if ind_lof.sum() > 10 else np.nan
    median_syn = df.loc[ind_syn, "raw_score"].median() if ind_syn.sum() > 10 else np.nan
    median_mis = df.loc[ind_mis, "raw_score"].median()

    temp = df["raw_score"].copy()

    if not np.isnan(median_syn) and not np.isnan(median_lof):
        if median_lof < median_syn:
            temp = -temp
            median_lof = -median_lof
            median_syn = -median_syn
        df["norm_raw_score"] = (temp - median_syn) / (median_lof - median_syn)

    elif not np.isnan(median_syn) and np.isnan(median_lof):
        mis_scores = df.loc[ind_mis, "raw_score"]
        if median_mis < median_syn:
            temp = -temp
            median_syn = -median_syn
            mis_scores = -mis_scores
        surrogate_lof = np.nanquantile(mis_scores, upper_bound)
        df["norm_raw_score"] = (temp - median_syn) / (surrogate_lof - median_syn)

    elif np.isnan(median_syn) and not np.isnan(median_lof):
        mis_scores = df.loc[ind_mis, "raw_score"]
        if median_mis > median_lof:
            temp = -temp
            median_lof = -median_lof
            mis_scores = -mis_scores
        surrogate_syn = np.nanquantile(mis_scores, lower_bound)
        df["norm_raw_score"] = (temp - surrogate_syn) / (median_lof - surrogate_syn)

    else:
        # No SYN or LOF: use Proline as orientation reference
        ind_pro = df["aaalt"] == "P"
        median_p = df.loc[ind_pro, "raw_score"].median()
        mis_scores = df.loc[ind_mis, "raw_score"]
        if median_p < median_mis:
            temp = -temp
            mis_scores = -mis_scores
        surrogate_syn = np.nanquantile(mis_scores, lower_bound)
        surrogate_lof = np.nanquantile(mis_scores, upper_bound)
        df["norm_raw_score"] = (temp - surrogate_syn) / (surrogate_lof - surrogate_syn)

    if force_quantile:
        norm_mis = df.loc[ind_mis, "norm_raw_score"]
        q_syn = np.nanquantile(norm_mis, lower_bound)
        q_lof = np.nanquantile(norm_mis, upper_bound)
        df["norm_raw_score"] = (df["norm_raw_score"] - q_syn) / (q_lof - q_syn)

    return df


# ---------------------------------------------------------------------------
# James-Stein shrinkage
# ---------------------------------------------------------------------------

def james_stein(matrix: np.ndarray) -> np.ndarray:
    """Compute James-Stein shrinkage estimates of column means.

    Mirrors ``get_js()`` from ``AA_proj_functions_101421.R``.

    Parameters
    ----------
    matrix:
        2-D array of shape (n_obs, n_positions).  NaN values are treated as
        missing (ignored in mean and variance calculations).

    Returns
    -------
    np.ndarray
        1-D array of length ``n_positions`` with shrinkage-adjusted means.
    """
    col_means = np.nanmean(matrix, axis=0)
    global_mean = np.nanmean(col_means)
    n_total = np.sum(~np.isnan(matrix))
    s2 = np.nanvar(matrix, ddof=0) / n_total if n_total > 0 else 0.0

    denom = np.sum((col_means - global_mean) ** 2)
    if denom == 0 or np.isnan(denom):
        return np.full(len(col_means), global_mean)

    p = len(col_means)
    cval = 1.0 - max(0.0, (p - 2) * s2 / denom)
    adj = cval * (col_means - global_mean)
    adj = np.where(np.isnan(adj), 0.0, adj)
    return global_mean + adj
