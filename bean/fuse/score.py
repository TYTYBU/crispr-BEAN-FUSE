"""
Core FUSE score computation.

Ports ``de_noise_ss_1gene()`` and the surrounding workflow from the original R
script (``LDLR_peDMS_FUSE_score.R`` / ``AA_proj_functions_101421.R``) to Python.

FUSE score = positional component + substitution component

* Positional component  – James-Stein shrinkage estimate of the mean effect at
  each amino-acid position across all observed substitutions.
* Substitution component – per-pair score from the FUNSUM substitution matrix.
  When a DSSP secondary-structure file is supplied, a secondary-structure-
  specific FUNSUM is used; otherwise the overall FUNSUM is applied.

The final score is ``pos_score + sub_score`` (or ``pos_score + sub_score_ss``
when secondary-structure information is available).
"""

from __future__ import annotations

import logging
import warnings
from typing import Optional

import numpy as np
import pandas as pd

from .dssp import parse_dssp
from .funsum import funsum_to_sub_table, load_funsum, load_funsum_ss
from .utils import (
    annotate_functional_class,
    extract_aa_columns,
    james_stein,
    normalize_scoreset,
    to_one_letter,
)

logger = logging.getLogger(__name__)

# DSSP codes → grouped secondary-structure label
_SS_GROUP: dict[str, str] = {
    "G": "Helices", "H": "Helices", "I": "Helices",
    "E": "Strands", "B": "Strands",
    "T": "Loops",   "S": "Loops",   "C": "Loops",
}


def _resolve_raw_score_col(df: pd.DataFrame, raw_score_col: Optional[str]) -> str:
    """Pick the raw-score column, preferring mu_z_adj, then mu_z."""
    if raw_score_col is not None:
        if raw_score_col not in df.columns:
            raise ValueError(
                f"Requested score column '{raw_score_col}' not found in input. "
                f"Available columns: {df.columns.tolist()}"
            )
        return raw_score_col

    if "mu_z_adj" in df.columns:
        logger.info("Using 'mu_z_adj' as raw score column.")
        return "mu_z_adj"

    if "mu_z" in df.columns:
        warnings.warn(
            "'mu_z_adj' not found in input. Falling back to 'mu_z' as the raw "
            "score.  Note: mu_z is not calibrated to a null distribution – "
            "interpret FUSE scores with caution.",
            UserWarning,
            stacklevel=3,
        )
        return "mu_z"

    raise ValueError(
        "Neither 'mu_z_adj' nor 'mu_z' found in input DataFrame. "
        "Please specify --score-col explicitly."
    )


def _build_position_matrix(
    df: pd.DataFrame,
    df_pos: pd.DataFrame,
    aa_list: list[str],
    df_funsum: pd.DataFrame,
    pos_mean_method: str,
) -> np.ndarray:
    """Build the (n_positions × n_aa) matrix of norm_raw_scores, then apply JS."""
    n_pos = len(df_pos)
    n_aa  = len(aa_list)
    aa_idx = {aa: i for i, aa in enumerate(aa_list)}

    mat = np.full((n_pos, n_aa), np.nan)

    # Optional residual matrix for pos_mean_method == "funsum"
    mat2 = np.full((n_pos, n_aa), np.nan) if pos_mean_method == "funsum" else None

    for i, row in df_pos.iterrows():
        gene, aapos = row["gene"], row["aapos"]
        sub = df[(df["gene"] == gene) & (df["aapos"] == aapos)]
        for _, srow in sub.iterrows():
            alt = srow["aaalt"]
            if alt in aa_idx:
                j = aa_idx[alt]
                mat[i, j] = srow["norm_raw_score"]
                if mat2 is not None:
                    ref = row["aaref"]
                    if ref in df_funsum.index and alt in df_funsum.columns:
                        mat2[i, j] = mat[i, j] - df_funsum.loc[ref, alt]

    return mat, mat2


def _compute_pos_means(
    mat: np.ndarray,
    mat2: Optional[np.ndarray],
    pos_mean_method: str,
) -> np.ndarray:
    """Compute per-position means using the chosen method.

    ``mat`` has shape (n_positions, n_aa).  Per-position means are computed
    across the AA axis (axis=1).  For James-Stein the matrix is transposed to
    (n_aa, n_positions) to match the R convention ``get_js(t(temp))``, where
    JS shrinkage is applied to the column means of that transposed matrix
    (= per-position means).
    """
    if pos_mean_method == "mean":
        return np.nanmean(mat, axis=1)
    elif pos_mean_method == "median":
        return np.nanmedian(mat, axis=1)
    elif pos_mean_method == "js":
        # R: get_js(t(temp))  →  colMeans of (n_aa × n_pos) = per-position means
        return james_stein(mat.T)
    elif pos_mean_method == "funsum":
        if mat2 is None:
            raise ValueError("mat2 required for pos_mean_method='funsum'")
        return james_stein(mat2.T)
    else:
        raise ValueError(f"Unknown pos_mean_method: '{pos_mean_method}'")


def calculate_fuse_scores(
    input_csv: str,
    edit_col: str = "target",
    raw_score_col: Optional[str] = None,
    mu_sd_max: Optional[float] = None,
    include_lof: bool = True,
    dss_path: Optional[str] = None,
    dssp_offset: int = 0,
    target_decoder: Optional[str] = None,
    pos_mean_method: str = "js",
    lower_bound: float = 0.1,
    upper_bound: float = 0.9,
    gene_name: Optional[str] = None,
) -> pd.DataFrame:
    """Compute FUSE scores from a bean element result CSV.

    Parameters
    ----------
    input_csv:
        Path to the ``bean_element_result.*.csv`` produced by ``bean run``.
    edit_col:
        Column containing amino-acid edit strings (e.g. ``"target"`` or
        ``"target_variant"``).  Strings are parsed as ``[gene:]pos:ref>alt``
        or the CodingNoncodingAllele ``aa_part|nt_part`` format.
    raw_score_col:
        Column to use as the raw score input to FUSE.  Defaults to
        ``mu_z_adj`` if present, otherwise ``mu_z`` (with a warning).
    mu_sd_max:
        If provided, rows where ``mu_sd > mu_sd_max`` are excluded before
        scoring.  Useful for filtering out low-confidence estimates.
    include_lof:
        If *True* (default), stop-gain (``*``) variants are included in the
        FUSE calculation using the LOF-inclusive FUNSUM.  If *False*, stop-gain
        variants are dropped and the noLOF FUNSUM is used.
    dss_path:
        Path to a DSSP ``.dss`` file for secondary-structure annotation.  When
        supplied, secondary-structure-specific FUNSUM scores are used to compute
        ``FUSE_SS_score``.  When absent, ``FUSE_score`` is also copied to
        ``FUSE_SS_score``.
    dssp_offset:
        Integer added to every PDB residue number (``respdb``) in the DSSP file
        before matching against ``aapos``.  Use when the PDB structure is
        numbered from the mature protein but ``aapos`` counts from the precursor
        (including signal peptide).  For example, ``dssp_offset=21`` shifts PDB
        residue 1 → ``aapos`` 22.  Default: ``0`` (no adjustment).
    target_decoder:
        Path to a CSV file that maps the ``target`` column to amino-acid
        coordinates.  Required columns: ``target``, ``aapos``, ``aaref``,
        ``aaalt`` (use ``Z`` or ``*`` for stop-gain).  Provide this when
        ``bean run`` was run without the ``bean filter`` AA-translation step,
        leaving ``target`` in a non-parseable format (e.g. ``16200_G_D``).
        When supplied, edit-string parsing is skipped entirely.
    pos_mean_method:
        Method for computing the positional component.  Fixed to ``"js"``
        (James-Stein shrinkage) by default.
    lower_bound, upper_bound:
        Quantile anchors used during score normalisation.
    gene_name:
        Optional gene label to attach to all variants (used in the
        ``gene_aa_str`` identifier column).  Inferred from the edit strings
        when not supplied.

    Returns
    -------
    pd.DataFrame
        All input rows (after filtering) with additional columns:

        ``gene``, ``aapos``, ``aaref``, ``aaalt``, ``functional_class``,
        ``raw_score``, ``norm_raw_score``,
        ``pos_score``, ``sub_score``, ``sub_score_ss`` (if DSSP available),
        ``FUSE_score``, ``FUSE_SS_score``, ``gene_aa_str``.
    """
    # ------------------------------------------------------------------
    # 1. Load data
    # ------------------------------------------------------------------
    df_raw = pd.read_csv(input_csv, index_col=0)

    # ------------------------------------------------------------------
    # 2. Resolve score column
    # ------------------------------------------------------------------
    score_col = _resolve_raw_score_col(df_raw, raw_score_col)
    logger.info(f"Using score column: '{score_col}'")

    # ------------------------------------------------------------------
    # 3. Optional mu_sd filter
    # ------------------------------------------------------------------
    if mu_sd_max is not None:
        if "mu_sd" not in df_raw.columns:
            warnings.warn(
                "'mu_sd' column not found – skipping mu_sd filter.",
                UserWarning,
                stacklevel=2,
            )
        else:
            n_before = len(df_raw)
            df_raw = df_raw[df_raw["mu_sd"] <= mu_sd_max].copy()
            logger.info(
                f"mu_sd filter (≤ {mu_sd_max}): kept {len(df_raw)}/{n_before} rows."
            )

    # ------------------------------------------------------------------
    # 4. Resolve aapos / aaref / aaalt
    #    Two paths:
    #      A. target_decoder CSV  – merge on "target" column (no string parsing)
    #      B. edit-string parsing – parse "[gene:]pos:ref>alt" from edit_col
    # ------------------------------------------------------------------
    df = df_raw.rename(columns={score_col: "raw_score"}).copy()

    if target_decoder is not None:
        logger.info(f"Loading target decoder from '{target_decoder}'.")
        df_dec = pd.read_csv(target_decoder)
        required_dec_cols = {"target", "aapos", "aaref", "aaalt"}
        missing_dec = required_dec_cols - set(df_dec.columns)
        if missing_dec:
            raise ValueError(
                f"target_decoder is missing required columns: {missing_dec}. "
                f"Found: {df_dec.columns.tolist()}"
            )
        # Normalise stop-codon encoding: Z → *
        def _norm_aa(x):
            try:
                return to_one_letter(str(x))
            except (ValueError, TypeError):
                return np.nan

        df_dec["aaref"] = df_dec["aaref"].map(_norm_aa)
        df_dec["aaalt"] = df_dec["aaalt"].map(_norm_aa)
        df_dec["aapos"] = pd.to_numeric(df_dec["aapos"], errors="coerce").astype("Int64")

        n_before = len(df)
        df = df.merge(
            df_dec[["target", "aapos", "aaref", "aaalt"]],
            on="target",
            how="left",
        )
        n_unmatched = df["aapos"].isna().sum()
        if n_unmatched:
            logger.info(
                f"{n_unmatched}/{n_before} rows had no match in target_decoder "
                "and will be excluded."
            )
    else:
        # Parse edit strings → aapos, aaref, aaalt
        df = extract_aa_columns(df, edit_col)

    # Drop rows with missing AA info or aapos ≤ 0
    valid = (
        df["aapos"].notna()
        & df["aaref"].notna()
        & df["aaalt"].notna()
        & (df["aapos"] > 0)
    )
    n_invalid = (~valid).sum()
    if n_invalid:
        logger.info(f"Dropping {n_invalid} rows with invalid / non-coding positions.")
    df = df[valid].copy()

    if df.empty:
        if target_decoder is not None:
            raise ValueError(
                "No valid amino-acid variants remain after merging target_decoder. "
                "Check that the 'target' values in the decoder match those in the input CSV."
            )
        raise ValueError(
            "No valid amino-acid variants remain after parsing the edit column. "
            f"Check that '{edit_col}' contains strings in the format "
            "'[gene:]pos:ref>alt' or 'aa_part|nt_part'.  "
            "If your target column uses a different format (e.g. '16200_G_D'), "
            "supply a target decoder CSV with --target-decoder."
        )

    # ------------------------------------------------------------------
    # 5. Infer gene name
    # ------------------------------------------------------------------
    if gene_name:
        df["gene"] = gene_name
    elif "gene" not in df.columns:
        if target_decoder is not None:
            df["gene"] = "gene"
        else:
            # Try to infer from edit strings like "LDLR:123:A>G"
            def _extract_gene(s: str) -> str:
                parts = s.split("|")[0].split(":")
                return parts[0] if len(parts) == 3 else "gene"

            df["gene"] = df[edit_col].map(_extract_gene)

    # ------------------------------------------------------------------
    # 6. Annotate functional class & filter LOF if needed
    # ------------------------------------------------------------------
    df = annotate_functional_class(df)

    if not include_lof:
        n_lof = (df["functional_class"] == "LOF").sum()
        if n_lof:
            logger.info(f"Dropping {n_lof} stop-gain (LOF) variants (include_lof=False).")
        df = df[~((df["aaref"] != "*") & (df["aaalt"] == "*"))].copy()
        df = df[df["aaref"] != "*"].copy()

    # ------------------------------------------------------------------
    # 7. Collapse duplicate substitutions (average raw_score)
    # ------------------------------------------------------------------
    df["gene_aa_str"] = df["gene"] + "---" + df["aaref"] + df["aapos"].astype(str) + df["aaalt"]
    df = (
        df.groupby("gene_aa_str", as_index=False)
        .agg(
            gene=("gene", "first"),
            aapos=("aapos", "first"),
            aaref=("aaref", "first"),
            aaalt=("aaalt", "first"),
            functional_class=("functional_class", "first"),
            raw_score=("raw_score", "mean"),
        )
    )

    # ------------------------------------------------------------------
    # 8. Normalise scores
    # ------------------------------------------------------------------
    df = normalize_scoreset(
        df,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        force_quantile=True,
    )

    # ------------------------------------------------------------------
    # 9. Load FUNSUM matrices
    # ------------------------------------------------------------------
    df_funsum    = load_funsum(include_lof=include_lof)
    ls_funsum_ss = load_funsum_ss(include_lof=include_lof)

    # Allowed alt AAs
    aa_list = list(df_funsum.columns)

    sub_table = funsum_to_sub_table(df_funsum)

    # ------------------------------------------------------------------
    # 10. Compute positional component
    # ------------------------------------------------------------------
    df_pos = (
        df.groupby(["gene", "aapos"], as_index=False)
        .agg(aaref=("aaref", "first"))
        .reset_index(drop=True)
    )
    df_pos["gene_aapos"] = df_pos["gene"] + "---" + df_pos["aapos"].astype(str)

    mat, mat2 = _build_position_matrix(df, df_pos, aa_list, df_funsum, pos_mean_method)
    df_pos["pos_mean"] = _compute_pos_means(mat, mat2, pos_mean_method)

    # ------------------------------------------------------------------
    # 11. Build full output table (all possible substitutions at observed pos)
    # ------------------------------------------------------------------
    rows = []
    for _, prow in df_pos.iterrows():
        for alt in aa_list:
            rows.append(
                {
                    "gene":       prow["gene"],
                    "aapos":      prow["aapos"],
                    "aaref":      prow["aaref"],
                    "aaalt":      alt,
                    "gene_aapos": prow["gene_aapos"],
                    "pos_score":  prow["pos_mean"],
                }
            )
    df_out = pd.DataFrame(rows)

    df_out["gene_aa_str"] = (
        df_out["gene"] + "---"
        + df_out["aaref"]
        + df_out["aapos"].astype(str)
        + df_out["aaalt"]
    )
    df_out["aa_pair"] = df_out["aaref"] + df_out["aaalt"]

    # Annotate functional class on the expanded table
    df_out = annotate_functional_class(df_out)

    # Merge in observed raw / normalised scores
    obs_cols = df[["gene_aa_str", "raw_score", "norm_raw_score"]]
    df_out = df_out.merge(obs_cols, on="gene_aa_str", how="left")

    # Substitution score (overall FUNSUM)
    sub_map = sub_table.set_index("aa_pair")["score"]
    df_out["sub_score"] = df_out["aa_pair"].map(sub_map)

    # ------------------------------------------------------------------
    # 12. Secondary-structure annotation (optional)
    # ------------------------------------------------------------------
    df_out["ss"]  = np.nan
    df_out["acc"] = np.nan
    df_out["sub_score_ss"] = np.nan

    if dss_path:
        try:
            df_dssp = parse_dssp(dss_path)
            # Apply offset: shift PDB residue numbers to match aapos numbering
            if dssp_offset != 0:
                logger.info(
                    f"Applying DSSP offset of {dssp_offset:+d} to all PDB residue "
                    f"numbers (respdb + {dssp_offset})."
                )
                df_dssp = df_dssp.copy()
                df_dssp["respdb"] = df_dssp["respdb"] + dssp_offset
            # Match by adjusted PDB residue number
            dssp_idx = df_dssp.set_index("respdb")
            df_out["ss"]  = df_out["aapos"].map(dssp_idx["ss"])
            df_out["acc"] = df_out["aapos"].map(dssp_idx["sasa"])

            for ss_code, ss_group in _SS_GROUP.items():
                mask = df_out["ss"] == ss_code
                if not mask.any():
                    continue
                if ss_group not in ls_funsum_ss:
                    continue
                ss_sub = funsum_to_sub_table(ls_funsum_ss[ss_group])
                ss_map = ss_sub.set_index("aa_pair")["score"]
                df_out.loc[mask, "sub_score_ss"] = df_out.loc[mask, "aa_pair"].map(ss_map)

        except Exception as exc:
            warnings.warn(
                f"Failed to parse DSSP file '{dss_path}': {exc}. "
                "Secondary-structure-specific scores will not be computed.",
                UserWarning,
                stacklevel=2,
            )

    # ------------------------------------------------------------------
    # 13. Compute final FUSE scores
    # ------------------------------------------------------------------
    df_out["FUSE_score"] = df_out["pos_score"] + df_out["sub_score"]

    # FUSE_SS_score: use SS-specific sub_score where available, fall back to overall
    df_out["FUSE_SS_score"] = df_out["pos_score"] + df_out["sub_score_ss"]
    missing_ss = df_out["FUSE_SS_score"].isna()
    df_out.loc[missing_ss, "FUSE_SS_score"] = df_out.loc[missing_ss, "FUSE_score"]
    df_out.loc[missing_ss, "sub_score_ss"]  = df_out.loc[missing_ss, "sub_score"]

    # ------------------------------------------------------------------
    # 14. Final column selection
    # ------------------------------------------------------------------
    col_order = [
        "gene", "aapos", "aaref", "aaalt", "functional_class",
        "raw_score", "norm_raw_score",
        "ss", "acc",
        "pos_score", "sub_score", "sub_score_ss",
        "FUSE_score", "FUSE_SS_score",
        "gene_aa_str",
    ]
    df_out = df_out[[c for c in col_order if c in df_out.columns]]

    return df_out
