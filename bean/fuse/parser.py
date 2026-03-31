"""
Argument parser for ``bean fuse``.
"""

from __future__ import annotations

import argparse


def add_fuse_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Populate *parser* with all ``bean fuse`` arguments and return it."""

    # ── Positional ────────────────────────────────────────────────────────────
    parser.add_argument(
        "input_csv",
        type=str,
        help=(
            "Path to a bean element result CSV produced by ``bean run`` "
            "(e.g. ``bean_element_result.*.csv``)."
        ),
    )

    # ── I/O options ───────────────────────────────────────────────────────────
    io_grp = parser.add_argument_group("Input / output options")

    io_grp.add_argument(
        "--edit-col",
        "-e",
        type=str,
        default="target",
        help=(
            "Column in the input CSV containing amino-acid edit strings in the "
            "format ``[gene:]pos:ref>alt`` or the CodingNoncodingAllele format "
            "``aa_part|nt_part``.  Ignored when ``--target-decoder`` is supplied.  "
            "Default: ``target``."
        ),
    )

    io_grp.add_argument(
        "--target-decoder",
        type=str,
        default=None,
        metavar="PATH",
        help=(
            "Path to a target-decoder CSV that maps the ``target`` column of the "
            "bean element result to amino-acid coordinates.  "
            "Required columns: ``target``, ``aapos`` (int), ``aaref`` (str), "
            "``aaalt`` (str, use ``Z`` or ``*`` for stop-gain).  "
            "Use this when ``bean run`` was executed without the ``bean filter`` "
            "AA-translation step, leaving the ``target`` column in a "
            "non-parseable format such as ``16200_G_D``."
        ),
    )

    io_grp.add_argument(
        "--score-col",
        "-s",
        type=str,
        default=None,
        help=(
            "Column to use as the raw score for FUSE.  "
            "Defaults to ``mu_z_adj`` if present in the CSV, otherwise falls "
            "back to ``mu_z`` with a warning."
        ),
    )

    io_grp.add_argument(
        "--outdir",
        "-o",
        type=str,
        default=".",
        help="Directory to write the output CSV.  Default: current directory.",
    )

    io_grp.add_argument(
        "--output-prefix",
        type=str,
        default=None,
        help=(
            "Prefix for the output file name.  "
            "Default: derived from the input file name."
        ),
    )

    # ── Filtering options ─────────────────────────────────────────────────────
    filt_grp = parser.add_argument_group("Filtering options")

    filt_grp.add_argument(
        "--mu-sd-max",
        type=float,
        default=None,
        metavar="FLOAT",
        help=(
            "Maximum allowed ``mu_sd`` value.  Variants with ``mu_sd`` above "
            "this threshold are excluded before FUSE scoring.  "
            "No filter is applied by default."
        ),
    )

    filt_grp.add_argument(
        "--include-lof",
        action="store_true",
        default=True,
        help=(
            "Include stop-gain (LoF / nonsense) variants in FUSE scoring using "
            "the LoF-inclusive FUNSUM matrix.  Enabled by default."
        ),
    )

    filt_grp.add_argument(
        "--exclude-lof",
        dest="include_lof",
        action="store_false",
        help=(
            "Exclude stop-gain (LoF / nonsense) variants from FUSE scoring and "
            "use the LoF-free FUNSUM matrix instead."
        ),
    )

    # ── Secondary structure ───────────────────────────────────────────────────
    ss_grp = parser.add_argument_group("Secondary structure options")

    ss_grp.add_argument(
        "--dss",
        type=str,
        default=None,
        metavar="PATH",
        help=(
            "Path to a DSSP output file (``.dss`` / ``.dssp``) generated from "
            "a PDB structure file using the DSSP program.  When provided, "
            "secondary-structure-specific FUNSUM scores are used to produce "
            "``FUSE_SS_score`` in addition to the overall ``FUSE_score``.  "
            "If omitted, ``FUSE_SS_score`` is set equal to ``FUSE_score``."
        ),
    )

    ss_grp.add_argument(
        "--dssp-offset",
        type=int,
        default=0,
        metavar="INT",
        help=(
            "Integer offset added to every PDB residue number (``respdb``) in "
            "the DSSP file before matching against the ``aapos`` values in the "
            "input CSV.  Use this when the PDB structure begins at the mature "
            "protein (after signal-peptide cleavage) but ``aapos`` is numbered "
            "from the full precursor sequence.  "
            "For example, if the LDLR signal peptide is 21 aa long and the PDB "
            "chain starts at residue 1 of the mature protein, pass "
            "``--dssp-offset 21`` so that PDB residue 1 maps to ``aapos`` 22.  "
            "Default: 0 (no adjustment)."
        ),
    )

    # ── Normalisation options ─────────────────────────────────────────────────
    norm_grp = parser.add_argument_group("Normalisation options")

    norm_grp.add_argument(
        "--lower-bound",
        type=float,
        default=0.1,
        metavar="FLOAT",
        help=(
            "Lower quantile of the missense distribution used as the SYN "
            "anchor when synonymous variants are absent.  Default: 0.1."
        ),
    )

    norm_grp.add_argument(
        "--upper-bound",
        type=float,
        default=0.9,
        metavar="FLOAT",
        help=(
            "Upper quantile of the missense distribution used as the LOF "
            "anchor when LoF variants are absent.  Default: 0.9."
        ),
    )

    norm_grp.add_argument(
        "--scale-fuse",
        dest="scale_fuse_scores",
        action="store_true",
        default=False,
        help=(
            "After computing FUSE_score and FUSE_SS_score, add normalised "
            "columns FUSE_score_norm and FUSE_SS_score_norm by rescaling so "
            "that the median of synonymous variants maps to 0 and the median "
            "of LoF variants maps to 1.  Requires both SYN and LOF variants "
            "to be present in the output; a warning is issued and the columns "
            "are omitted if either group is absent."
        ),
    )

    # ── Misc ──────────────────────────────────────────────────────────────────
    misc_grp = parser.add_argument_group("Miscellaneous")

    misc_grp.add_argument(
        "--gene",
        type=str,
        default=None,
        metavar="STR",
        help=(
            "Gene label to attach to all variants.  Inferred from the edit "
            "strings when not supplied (e.g. ``LDLR`` from ``LDLR:123:A>G``)."
        ),
    )

    misc_grp.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        default=False,
        help="Print detailed progress messages.",
    )

    return parser


def get_fuse_parser(
    parent: argparse.ArgumentParser | None = None,
) -> argparse.ArgumentParser:
    """Return a fully configured argument parser for ``bean fuse``."""
    if parent is None:
        parent = argparse.ArgumentParser(
            prog="bean fuse",
            description=(
                "Compute FUSE (Functional Score Using Structural Ensemble) "
                "scores from a bean element result CSV.\n\n"
                "FUSE denoises per-variant functional scores by combining a "
                "James-Stein positional estimate with a substitution-specific "
                "score from the FUNSUM substitution matrix.  Optionally, "
                "secondary-structure-specific FUNSUM matrices are applied when "
                "a DSSP file is provided."
            ),
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
    return add_fuse_args(parent)
