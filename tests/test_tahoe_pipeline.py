import unittest

import numpy as np

from src.tahoe_pipeline import build_sparse_expression_matrix, normalize_cpm_log1p


class TahoePipelineTests(unittest.TestCase):
    def test_build_sparse_expression_matrix_drops_sentinel_and_preserves_metadata(self):
        stream = [
            {
                "genes": [-1, 101, 202, 999],
                "expressions": [-1.0, 3.0, 7.0, 11.0],
                "drug": "compound-a",
            },
            {
                "genes": [101],
                "expressions": [5.0],
                "drug": "compound-b",
            },
        ]
        gene_vocab = {101: "ENSG_A", 202: "ENSG_B"}

        matrix, obs_df, var_df = build_sparse_expression_matrix(stream, gene_vocab, n_cells=2)

        self.assertEqual(matrix.shape, (2, 2))
        self.assertEqual(matrix.nnz, 3)
        self.assertEqual(obs_df["drug"].tolist(), ["compound-a", "compound-b"])
        self.assertEqual(var_df["ensembl_id"].tolist(), ["ENSG_A", "ENSG_B"])
        np.testing.assert_array_equal(matrix.toarray(), np.array([[3.0, 7.0], [5.0, 0.0]]))

    def test_normalize_cpm_log1p_handles_zero_rows(self):
        stream = [
            {"genes": [101, 202], "expressions": [1.0, 3.0]},
            {"genes": [], "expressions": []},
        ]
        gene_vocab = {101: "ENSG_A", 202: "ENSG_B"}
        matrix, _, _ = build_sparse_expression_matrix(stream, gene_vocab, n_cells=2)

        normalized = normalize_cpm_log1p(matrix).toarray()

        expected_first_row = np.log1p(np.array([250_000.0, 750_000.0]))
        np.testing.assert_allclose(normalized[0], expected_first_row)
        np.testing.assert_array_equal(normalized[1], np.array([0.0, 0.0]))


if __name__ == "__main__":
    unittest.main()
