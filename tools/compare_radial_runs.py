#!/usr/bin/env python3
"""Compare two archived radial solver runs on the reference radial grid."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np


class ComparisonError(RuntimeError):
    pass


def read_table(path: Path, minimum_columns: int) -> np.ndarray:
    try:
        data = np.loadtxt(path, comments="#", ndmin=2)
    except (OSError, ValueError) as exc:
        raise ComparisonError(f"could not read numeric data from {path}") from exc
    if data.shape[1] < minimum_columns:
        raise ComparisonError(
            f"{path} has {data.shape[1]} columns; expected at least {minimum_columns}"
        )
    if not np.all(np.isfinite(data)):
        raise ComparisonError(f"{path} contains a non-finite value")
    if np.any(np.diff(data[:, 0]) <= 0.0):
        raise ComparisonError(f"the radial coordinate is not strictly increasing in {path}")
    return data


def interpolate_to_reference(
    reference: np.ndarray, candidate: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    radius = reference[:, 0]
    tolerance = 10.0 * np.finfo(float).eps * max(1.0, float(np.max(np.abs(radius))))
    if radius[0] < candidate[0, 0] - tolerance or radius[-1] > candidate[-1, 0] + tolerance:
        raise ComparisonError(
            "candidate radial domain does not cover the reference radial domain"
        )
    values = np.column_stack(
        [np.interp(radius, candidate[:, 0], candidate[:, column])
         for column in range(1, candidate.shape[1])]
    )
    return reference[:, 1:], values


def compare_values(
    reference: np.ndarray,
    candidate: np.ndarray,
    absolute_tolerance: float,
    relative_tolerance: float,
) -> Dict[str, object]:
    if reference.shape != candidate.shape:
        raise ComparisonError(
            f"comparison shapes differ: {reference.shape} and {candidate.shape}"
        )
    difference = candidate - reference
    scale = absolute_tolerance + relative_tolerance * np.abs(reference)
    scaled_error = np.abs(difference) / scale
    reference_norm = float(np.linalg.norm(reference.ravel()))
    difference_norm = float(np.linalg.norm(difference.ravel()))
    return {
        "passed": bool(np.all(scaled_error <= 1.0)),
        "maximum_absolute_error": float(np.max(np.abs(difference))),
        "rms_absolute_error": float(np.sqrt(np.mean(np.abs(difference) ** 2))),
        "relative_l2_error": (
            difference_norm / reference_norm if reference_norm > 0.0 else None
        ),
        "maximum_tolerance_scaled_error": float(np.max(scaled_error)),
        "values_compared": int(reference.size),
    }


def normalized_current_columns(data: np.ndarray) -> np.ndarray:
    if data.shape[1] >= 6:
        # Legacy: r, pair amplitude, transverse magnitude, vx, vy, vz.
        return np.column_stack((data[:, 0], data[:, 1], data[:, 3:6]))
    # new_src: r, pair amplitude, vx, vy, vz.
    return data[:, :5]


def compare_file(
    reference_path: Path,
    candidate_path: Path,
    minimum_columns: int,
    absolute_tolerance: float,
    relative_tolerance: float,
    normalize_current: bool = False,
) -> Dict[str, object]:
    reference = read_table(reference_path, minimum_columns)
    candidate = read_table(candidate_path, minimum_columns)
    if normalize_current:
        reference = normalized_current_columns(reference)
        candidate = normalized_current_columns(candidate)
    reference_values, candidate_values = interpolate_to_reference(reference, candidate)
    result = compare_values(
        reference_values, candidate_values, absolute_tolerance, relative_tolerance
    )
    result.update(
        {
            "reference_points": int(reference.shape[0]),
            "candidate_points": int(candidate.shape[0]),
            "radial_interval": [float(reference[0, 0]), float(reference[-1, 0])],
            "candidate_interpolated_to_reference_grid": not (
                reference.shape[0] == candidate.shape[0]
                and np.array_equal(reference[:, 0], candidate[:, 0])
            ),
        }
    )
    return result


def parse_arguments(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare op_xyz, op_harm, and curr in two radial run directories. "
            "The candidate is interpolated component-wise onto the reference grid."
        )
    )
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--atol", type=float, default=1.0e-10)
    parser.add_argument("--rtol", type=float, default=1.0e-8)
    parser.add_argument("--output", type=Path, help="optional JSON report path")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_arguments(argv)
    if args.atol < 0.0 or args.rtol < 0.0:
        print("error: tolerances must be non-negative", file=sys.stderr)
        return 2
    reference_dir = args.reference.expanduser().resolve()
    candidate_dir = args.candidate.expanduser().resolve()
    try:
        comparisons = {
            "op_xyz": compare_file(
                reference_dir / "op_xyz", candidate_dir / "op_xyz", 19,
                args.atol, args.rtol
            ),
            "op_harm": compare_file(
                reference_dir / "op_harm", candidate_dir / "op_harm", 19,
                args.atol, args.rtol
            ),
            "curr": compare_file(
                reference_dir / "curr", candidate_dir / "curr", 5,
                args.atol, args.rtol, normalize_current=True
            ),
        }
    except ComparisonError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    passed = all(bool(result["passed"]) for result in comparisons.values())
    report: Dict[str, object] = {
        "passed": passed,
        "reference_directory": str(reference_dir),
        "candidate_directory": str(candidate_dir),
        "absolute_tolerance": args.atol,
        "relative_tolerance": args.rtol,
        "comparison_grid": "reference",
        "gauge_alignment": "none; inputs must use the same gauge and orientation",
        "comparisons": comparisons,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        output_path = args.output.expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
