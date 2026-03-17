#!/usr/bin/env python3

from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from validation.run_validation_suite import read_vtk_velocity, relative_l2_difference


class ValidationSuiteHelperTests(unittest.TestCase):
    def test_relative_l2_difference_matches_expected_value(self) -> None:
        value = relative_l2_difference([2.0, 4.0], [1.0, 2.0])
        self.assertTrue(math.isclose(value, 1.0, rel_tol=0.0, abs_tol=1.0e-12))

    def test_relative_l2_difference_rejects_mismatched_lengths(self) -> None:
        with self.assertRaisesRegex(
            ValueError, "equal-length vectors, got 2 and 1"
        ):
            relative_l2_difference([1.0, 1000.0], [1.0])

    def test_relative_l2_difference_rejects_empty_vectors(self) -> None:
        with self.assertRaisesRegex(ValueError, "non-empty vectors"):
            relative_l2_difference([], [])

    def test_read_vtk_velocity_parses_expected_vector_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            vtk_path = Path(tmpdir) / "valid.vtk"
            vtk_path.write_text(
                "\n".join(
                    [
                        "# vtk DataFile Version 3.0",
                        "test",
                        "ASCII",
                        "DATASET STRUCTURED_POINTS",
                        "DIMENSIONS 1 1 1",
                        "POINT_DATA 1",
                        "VECTORS velocity double",
                        "1 2 3",
                        "SCALARS pressure double 1",
                        "LOOKUP_TABLE default",
                        "0",
                    ]
                )
                + "\n"
            )

            values = read_vtk_velocity(vtk_path, 3)
            self.assertEqual(values, [1.0, 2.0, 3.0])

    def test_read_vtk_velocity_rejects_truncated_field(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            vtk_path = Path(tmpdir) / "truncated.vtk"
            vtk_path.write_text(
                "\n".join(
                    [
                        "# vtk DataFile Version 3.0",
                        "test",
                        "ASCII",
                        "DATASET STRUCTURED_POINTS",
                        "DIMENSIONS 1 1 1",
                        "POINT_DATA 1",
                        "VECTORS velocity double",
                        "1 2 3",
                        "SCALARS pressure double 1",
                        "LOOKUP_TABLE default",
                        "0",
                    ]
                )
                + "\n"
            )

            with self.assertRaisesRegex(ValueError, "expected 6 values"):
                read_vtk_velocity(vtk_path, 6)


if __name__ == "__main__":
    unittest.main()
