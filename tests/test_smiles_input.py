import unittest
from pathlib import Path

from conformer_app.core.config import EmbedConfig, MMConfig, OutputConfig
from conformer_app.core.runner import run_mm_for_job


class TestSmilesInput(unittest.TestCase):
    def setUp(self) -> None:
        self.output_dir = Path("sample") / "test1"
        if self.output_dir.exists():
            for entry in self.output_dir.iterdir():
                if entry.name == ".gitkeep":
                    continue
                if entry.is_dir():
                    for sub in entry.rglob("*"):
                        if sub.is_file():
                            sub.unlink()
                    for sub in sorted(entry.rglob("*"), reverse=True):
                        if sub.is_dir():
                            sub.rmdir()
                    entry.rmdir()
                else:
                    entry.unlink()
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def test_smiles_input_generates_post_mm_xyz(self) -> None:
        embed_cfg = EmbedConfig(num_confs=5)
        mm_cfg = MMConfig(max_keep=2)
        out_cfg = OutputConfig(base_dir=str(self.output_dir))

        result = run_mm_for_job("mol1", "C1=CC=CC=C1", embed_cfg, mm_cfg, out_cfg)

        self.assertEqual(result.n_selected, 2)
        for xyz_path in result.xyz_files:
            self.assertTrue(xyz_path.exists())
            content = xyz_path.read_text(encoding="utf-8")
            self.assertIn("C1=CC=CC=C1", content)


if __name__ == "__main__":
    unittest.main()
