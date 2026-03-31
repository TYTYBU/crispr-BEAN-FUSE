"""
CLI entry point for ``bean fuse``.

Computes FUSE (Functional Score Using Structural Ensemble) scores from a
bean element result CSV and writes the output to a CSV file.
"""

from __future__ import annotations

import logging
import os
import sys


def main(args) -> None:
    """Run FUSE scoring from parsed CLI arguments."""
    # ── Logging setup ─────────────────────────────────────────────────────────
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [bean fuse] %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
        level=level,
        stream=sys.stderr,
    )
    logger = logging.getLogger(__name__)

    # ── Validate inputs ────────────────────────────────────────────────────────
    if not os.path.isfile(args.input_csv):
        logger.error(f"Input file not found: {args.input_csv}")
        sys.exit(1)

    if args.dss and not os.path.isfile(args.dss):
        logger.error(f"DSSP file not found: {args.dss}")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    # ── Determine output path ──────────────────────────────────────────────────
    if args.output_prefix:
        prefix = args.output_prefix
    else:
        base = os.path.basename(args.input_csv)
        # e.g. bean_element_result.MixtureNormal.csv → bean_element_result.MixtureNormal
        prefix = os.path.splitext(base)[0]

    out_path = os.path.join(args.outdir, f"{prefix}.FUSE_score.csv")

    # ── Run FUSE scoring ───────────────────────────────────────────────────────
    logger.info(f"Input:         {args.input_csv}")
    logger.info(f"Edit column:   {args.edit_col}")
    if args.score_col:
        logger.info(f"Score column:  {args.score_col}")
    if args.mu_sd_max is not None:
        logger.info(f"mu_sd filter:  ≤ {args.mu_sd_max}")
    logger.info(f"Include LoF:   {args.include_lof}")
    if args.target_decoder:
        logger.info(f"Target decoder:{args.target_decoder}")
    if args.dss:
        logger.info(f"DSSP file:     {args.dss}")
        if args.dssp_offset:
            logger.info(f"DSSP offset:   {args.dssp_offset:+d}")
    else:
        logger.info("DSSP file:     not provided (using overall FUNSUM only)")
    if args.scale_fuse_scores:
        logger.info("Scale FUSE:    enabled (SYN=0, LOF=1 normalisation)")
    logger.info(f"Output:        {out_path}")

    from bean.fuse.score import calculate_fuse_scores

    try:
        df_fuse = calculate_fuse_scores(
            input_csv=args.input_csv,
            edit_col=args.edit_col,
            raw_score_col=args.score_col,
            mu_sd_max=args.mu_sd_max,
            include_lof=args.include_lof,
            dss_path=args.dss,
            dssp_offset=args.dssp_offset,
            target_decoder=args.target_decoder,
            lower_bound=args.lower_bound,
            upper_bound=args.upper_bound,
            gene_name=args.gene,
            scale_fuse_scores=args.scale_fuse_scores,
        )
    except Exception as exc:
        logger.error(f"FUSE scoring failed: {exc}")
        raise

    df_fuse.to_csv(out_path, index=False)
    logger.info(
        f"Done. {len(df_fuse)} variant×substitution records written to {out_path}"
    )
