#!/usr/bin/env python3
"""Regression check for the BH cold-channel AM suppression convention."""

from __future__ import annotations

import math
import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
SOURCE = REPO_ROOT / "src/enzo/Grid_BHAccretionHandler.C"


def expected_f_am(cs: float, vrot: float, cvisc: float) -> float:
    if vrot <= 0.0:
        return 1.0
    return min(1.0, (cs / vrot) ** 3 / cvisc)


def source_f_am_expression() -> str:
    text = SOURCE.read_text(encoding="utf-8")
    match = re.search(
        r"f_am\s*=\s*min\s*\(\s*1\.0\s*,\s*(.*?)\s*\);",
        text,
        re.DOTALL,
    )
    if match is None:
        raise AssertionError("Could not find f_am min(...) assignment in source")
    return " ".join(match.group(1).split())


def test_expected_reference_values() -> None:
    cvisc = 2.0 * math.pi
    assert math.isclose(expected_f_am(1.0, 1.0, cvisc), 1.0 / cvisc)
    assert math.isclose(expected_f_am(1.0, 2.0, cvisc), 0.125 / cvisc)
    assert expected_f_am(10.0, 1.0, cvisc) == 1.0
    assert expected_f_am(1.0, 0.0, cvisc) == 1.0


def test_source_uses_inverse_cvisc() -> None:
    expression = source_f_am_expression()
    if "BHAccretionCVisc) * pow" in expression:
        raise AssertionError(f"f_am multiplies by CVisc instead of dividing: {expression}")
    expected = "pow(cold_cs_avg / v_rot_cold, 3.0) / double(BHAccretionCVisc)"
    assert expected in expression, expression


def main() -> None:
    test_expected_reference_values()
    test_source_uses_inverse_cvisc()
    print("BH accretion AM suppression formula checks passed")


if __name__ == "__main__":
    main()
