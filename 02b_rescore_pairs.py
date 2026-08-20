#!/usr/bin/env python3
"""
03_rescore_pairs.py

Re-score primer pairs that have already passed the hard filters.

Input:
    pairs_pass.tsv
        Produced by 02_pair_primers.py.

Output:
    pairs_rescored.tsv
        All primer pairs sorted from best to worst by score_preblast.
        Lower score_preblast is better.

Scoring principle:
    1. Hard filters are handled upstream.
    2. Every soft criterion is converted to a penalty in [0, 1].
       0 = ideal / safe
       1 = close to the allowed limit / least preferred
    3. Weighted penalties are summed.
    4. Genome specificity is NOT included here. BLAST specificity should
       remain a downstream hard biological criterion.

This script uses only the Python standard library.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List


# ============================================================
# 1. Scoring parameters
# ============================================================

# Primer Tm
TM_MIN = 57.0
TM_OPT = 60.0
TM_MAX = 63.0

# Primer GC%
GC_MIN = 40.0
GC_OPT = 50.0
GC_MAX = 60.0

# Maximum allowed F/R Tm difference upstream
DELTA_TM_MAX = 2.5

# Secondary-structure penalty:
# <= STRUCTURE_SAFE: no penalty
# >= STRUCTURE_CUTOFF: penalty = 1
STRUCTURE_SAFE = 35.0
STRUCTURE_CUTOFF = 47.0

# Product size.
# The pair-generation hard range remains 50-700 bp.
# For ordinary endpoint PCR / gel scoring, 150-400 bp is treated as
# an equally good plateau rather than forcing a unique optimum of 250 bp.
PRODUCT_MIN = 50
PRODUCT_PREFERRED_MIN = 150
PRODUCT_PREFERRED_MAX = 400
PRODUCT_MAX = 700


# ============================================================
# 2. Weights
# ============================================================
#
# All component penalties are normalized to [0, 1].
# Weights sum to 1.0.
#
# Priority for a conventional sdY endpoint-PCR marker:
# pair 3'-interaction > heterodimer > delta-Tm > Tm >
# self 3'-interaction > hairpin > homodimer >
# 3'-GC > whole-primer GC > product size
#
WEIGHTS = {
    "pair_end": 0.20,
    "heterodimer": 0.16,
    "delta_tm": 0.15,
    "tm": 0.12,
    "self_end": 0.10,
    "hairpin": 0.08,
    "homodimer": 0.07,
    "end_gc": 0.05,
    "gc": 0.04,
    "product": 0.03,
}


# Output columns added by this script.
SCORE_FIELDS = [
    "penalty_pair_end",
    "penalty_heterodimer",
    "penalty_delta_tm",
    "penalty_tm",
    "penalty_self_end",
    "penalty_hairpin",
    "penalty_homodimer",
    "penalty_end_gc",
    "penalty_gc",
    "penalty_product",
    "score_preblast",
    "score_preblast_100",
    "rank_preblast",
]


REQUIRED_FIELDS = {
    "F_tm",
    "R_tm",
    "delta_tm",
    "F_gc",
    "R_gc",
    "F_end_gc5",
    "R_end_gc5",
    "F_hairpin_tm",
    "R_hairpin_tm",
    "F_homodimer_tm",
    "R_homodimer_tm",
    "F_self_end_tm",
    "R_self_end_tm",
    "heterodimer_tm",
    "pair_end_tm",
    "product_size",
}


# ============================================================
# 3. Utility functions
# ============================================================

def clamp01(x: float) -> float:
    """Clamp a number to the interval [0, 1]."""
    return max(0.0, min(1.0, x))


def require_finite(value: float, field_name: str) -> float:
    """Reject NaN/inf values explicitly."""
    if not math.isfinite(value):
        raise ValueError(f"{field_name} is not finite: {value}")
    return value


def as_float(row: Dict[str, str], key: str) -> float:
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid numeric value for column '{key}': {row.get(key)!r}"
        ) from exc
    return require_finite(value, key)


def as_int(row: Dict[str, str], key: str) -> int:
    try:
        # float() also accepts values such as "250.0".
        value = float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid integer-like value for column '{key}': {row.get(key)!r}"
        ) from exc

    require_finite(value, key)

    if not value.is_integer():
        raise ValueError(
            f"Column '{key}' should be an integer-like value, got {value}"
        )
    return int(value)


def two_sided_penalty(
    value: float,
    optimum: float,
    lower_limit: float,
    upper_limit: float,
) -> float:
    """
    Convert a value to [0,1] penalty around an optimum.

    optimum -> 0
    lower_limit or upper_limit -> 1

    Values beyond the limits are clamped to 1. This script normally receives
    only pairs that already passed hard filters, so beyond-limit values should
    be unusual.
    """
    if value == optimum:
        return 0.0

    if value < optimum:
        denom = optimum - lower_limit
        if denom <= 0:
            raise ValueError("Invalid lower scoring range.")
        return clamp01((optimum - value) / denom)

    denom = upper_limit - optimum
    if denom <= 0:
        raise ValueError("Invalid upper scoring range.")
    return clamp01((value - optimum) / denom)


def structure_penalty(
    tm: float,
    safe: float = STRUCTURE_SAFE,
    cutoff: float = STRUCTURE_CUTOFF,
) -> float:
    """
    Penalize predicted secondary-structure Tm.

    tm <= safe   -> 0
    safe < tm < cutoff -> linear 0..1
    tm >= cutoff -> 1
    """
    if cutoff <= safe:
        raise ValueError("STRUCTURE_CUTOFF must be greater than STRUCTURE_SAFE.")

    if tm <= safe:
        return 0.0

    return clamp01((tm - safe) / (cutoff - safe))


def end_gc_penalty(gc_count_last5: int) -> float:
    """
    Soft preference for a moderate GC content in the last 5 nt.

    last-5 GC count:
        1, 2, 3 -> 0.0
        0, 4    -> 0.5
        5       -> 1.0

    This is intentionally a broad plateau, not a rule that exactly
    two GC bases must occur at the 3' end.
    """
    if not 0 <= gc_count_last5 <= 5:
        raise ValueError(
            f"GC count in last 5 nt must be 0..5; got {gc_count_last5}"
        )

    if 1 <= gc_count_last5 <= 3:
        return 0.0

    if gc_count_last5 in (0, 4):
        return 0.5

    return 1.0


def product_penalty(size: int) -> float:
    """
    Product-size preference for conventional PCR + agarose gel.

    150-400 bp -> penalty 0

    50-150 bp:
        penalty increases linearly from 0 at 150 to 1 at 50

    400-700 bp:
        penalty increases linearly from 0 at 400 to 1 at 700
    """
    if PRODUCT_PREFERRED_MIN <= size <= PRODUCT_PREFERRED_MAX:
        return 0.0

    if size < PRODUCT_PREFERRED_MIN:
        denom = PRODUCT_PREFERRED_MIN - PRODUCT_MIN
        if denom <= 0:
            raise ValueError("Invalid lower product-size scoring range.")
        return clamp01(
            (PRODUCT_PREFERRED_MIN - size) / denom
        )

    denom = PRODUCT_MAX - PRODUCT_PREFERRED_MAX
    if denom <= 0:
        raise ValueError("Invalid upper product-size scoring range.")

    return clamp01(
        (size - PRODUCT_PREFERRED_MAX) / denom
    )


# ============================================================
# 4. Score one primer pair
# ============================================================

def score_pair(row: Dict[str, str]) -> Dict[str, float]:
    """
    Calculate all normalized component penalties and final weighted score.

    Returns a dict whose penalties are all in [0,1].
    score_preblast is also in [0,1] because weights sum to 1.
    """
    f_tm = as_float(row, "F_tm")
    r_tm = as_float(row, "R_tm")

    # Recalculate delta-Tm rather than trusting the stored rounded value.
    delta_tm = abs(f_tm - r_tm)

    f_gc = as_float(row, "F_gc")
    r_gc = as_float(row, "R_gc")

    f_end_gc5 = as_int(row, "F_end_gc5")
    r_end_gc5 = as_int(row, "R_end_gc5")

    f_hairpin_tm = as_float(row, "F_hairpin_tm")
    r_hairpin_tm = as_float(row, "R_hairpin_tm")

    f_homodimer_tm = as_float(row, "F_homodimer_tm")
    r_homodimer_tm = as_float(row, "R_homodimer_tm")

    f_self_end_tm = as_float(row, "F_self_end_tm")
    r_self_end_tm = as_float(row, "R_self_end_tm")

    heterodimer_tm = as_float(row, "heterodimer_tm")
    pair_end_tm = as_float(row, "pair_end_tm")

    product_size = as_int(row, "product_size")

    # ----- Tm -----
    p_tm_f = two_sided_penalty(
        f_tm, TM_OPT, TM_MIN, TM_MAX
    )
    p_tm_r = two_sided_penalty(
        r_tm, TM_OPT, TM_MIN, TM_MAX
    )
    p_tm = (p_tm_f + p_tm_r) / 2.0

    # ----- delta Tm -----
    p_delta_tm = clamp01(delta_tm / DELTA_TM_MAX)

    # ----- whole-primer GC -----
    p_gc_f = two_sided_penalty(
        f_gc, GC_OPT, GC_MIN, GC_MAX
    )
    p_gc_r = two_sided_penalty(
        r_gc, GC_OPT, GC_MIN, GC_MAX
    )
    p_gc = (p_gc_f + p_gc_r) / 2.0

    # ----- 3' GC -----
    p_end_gc = (
        end_gc_penalty(f_end_gc5)
        + end_gc_penalty(r_end_gc5)
    ) / 2.0

    # ----- hairpins -----
    p_hairpin = (
        structure_penalty(f_hairpin_tm)
        + structure_penalty(r_hairpin_tm)
    ) / 2.0

    # ----- homodimers -----
    p_homodimer = (
        structure_penalty(f_homodimer_tm)
        + structure_penalty(r_homodimer_tm)
    ) / 2.0

    # ----- self 3'-end interaction -----
    p_self_end = (
        structure_penalty(f_self_end_tm)
        + structure_penalty(r_self_end_tm)
    ) / 2.0

    # ----- F/R heterodimer -----
    p_heterodimer = structure_penalty(heterodimer_tm)

    # ----- F/R 3'-end interaction -----
    p_pair_end = structure_penalty(pair_end_tm)

    # ----- product size -----
    p_product = product_penalty(product_size)

    penalties = {
        "pair_end": p_pair_end,
        "heterodimer": p_heterodimer,
        "delta_tm": p_delta_tm,
        "tm": p_tm,
        "self_end": p_self_end,
        "hairpin": p_hairpin,
        "homodimer": p_homodimer,
        "end_gc": p_end_gc,
        "gc": p_gc,
        "product": p_product,
    }

    score = sum(
        WEIGHTS[name] * penalties[name]
        for name in WEIGHTS
    )

    return {
        "penalty_pair_end": p_pair_end,
        "penalty_heterodimer": p_heterodimer,
        "penalty_delta_tm": p_delta_tm,
        "penalty_tm": p_tm,
        "penalty_self_end": p_self_end,
        "penalty_hairpin": p_hairpin,
        "penalty_homodimer": p_homodimer,
        "penalty_end_gc": p_end_gc,
        "penalty_gc": p_gc,
        "penalty_product": p_product,
        "score_preblast": score,
        "score_preblast_100": 100.0 * score,
    }


# ============================================================
# 5. I/O
# ============================================================

def validate_weights() -> None:
    total = sum(WEIGHTS.values())

    if any(w < 0 for w in WEIGHTS.values()):
        raise ValueError("All weights must be non-negative.")

    if not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError(
            f"WEIGHTS must sum to 1.0; currently sum to {total:.12f}"
        )


def validate_header(fieldnames: List[str] | None) -> None:
    if not fieldnames:
        raise ValueError("Input TSV has no header.")

    missing = sorted(REQUIRED_FIELDS - set(fieldnames))
    if missing:
        raise ValueError(
            "Input TSV is missing required columns:\n  "
            + "\n  ".join(missing)
        )


def load_and_score(path: str | Path) -> tuple[List[Dict[str, str]], List[str]]:
    rows: List[Dict[str, str]] = []

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        validate_header(reader.fieldnames)

        original_fields = list(reader.fieldnames or [])

        for line_no, row in enumerate(reader, start=2):
            try:
                scores = score_pair(row)
            except Exception as exc:
                pair_id = row.get("pair_id", "NA")
                raise ValueError(
                    f"Failed to score line {line_no}, pair {pair_id}: {exc}"
                ) from exc

            # Preserve the previous score if this file has already been scored.
            if "score_preblast" in row and "score_preblast_old" not in row:
                row["score_preblast_old"] = row["score_preblast"]

            for key, value in scores.items():
                row[key] = f"{value:.6f}"

            rows.append(row)

    # Lower score is better.
    # pair_id is used only as a deterministic tie-breaker.
    rows.sort(
        key=lambda r: (
            float(r["score_preblast"]),
            r.get("pair_id", ""),
        )
    )

    for rank, row in enumerate(rows, start=1):
        row["rank_preblast"] = str(rank)

    # Preserve all original columns, add score_preblast_old if needed,
    # then append the new score components.
    output_fields: List[str] = []

    for field in original_fields:
        # score_preblast will be appended with the new score fields below.
        if field != "score_preblast":
            output_fields.append(field)

    if "score_preblast" in original_fields:
        output_fields.append("score_preblast_old")

    output_fields.extend(SCORE_FIELDS)

    # Remove accidental duplicates while preserving order.
    output_fields = list(dict.fromkeys(output_fields))

    return rows, output_fields


def write_tsv(
    rows: Iterable[Dict[str, str]],
    fieldnames: List[str],
    output_path: str | Path,
) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def print_summary(rows: List[Dict[str, str]], top_n: int = 10) -> None:
    print(f"Primer pairs scored: {len(rows):,}")

    if not rows:
        return

    scores = [float(row["score_preblast"]) for row in rows]

    print(f"Best score:   {min(scores):.6f}")
    print(f"Worst score:  {max(scores):.6f}")
    print()

    print(f"Top {min(top_n, len(rows))} primer pairs:")
    print(
        "rank\tpair_id\tscore\tF_tm\tR_tm\tdelta_tm\t"
        "product_size\tpair_end_penalty\theterodimer_penalty"
    )

    for row in rows[:top_n]:
        print(
            f"{row['rank_preblast']}\t"
            f"{row.get('pair_id', 'NA')}\t"
            f"{row['score_preblast']}\t"
            f"{row.get('F_tm', 'NA')}\t"
            f"{row.get('R_tm', 'NA')}\t"
            f"{row.get('delta_tm', 'NA')}\t"
            f"{row.get('product_size', 'NA')}\t"
            f"{row['penalty_pair_end']}\t"
            f"{row['penalty_heterodimer']}"
        )


# ============================================================
# 6. Main
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Re-score PASS primer pairs with normalized 0-1 penalties. "
            "Lower score_preblast is better."
        )
    )

    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Input pairs_pass.tsv from 02_pair_primers.py",
    )

    parser.add_argument(
        "--output",
        "-o",
        default="out/pairs_rescored.tsv",
        help="Output TSV sorted by score_preblast (default: out/pairs_rescored.tsv)",
    )

    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Number of top-ranked pairs to print to stdout (default: 10)",
    )

    args = parser.parse_args()

    validate_weights()

    rows, fieldnames = load_and_score(args.input)
    write_tsv(rows, fieldnames, args.output)

    print_summary(rows, top_n=max(0, args.top))
    print()
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
