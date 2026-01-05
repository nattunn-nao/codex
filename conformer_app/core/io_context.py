from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import pandas as pd

from .status import JobProgressRow


@dataclass
class IOContext:
    """
    ファイル入出力に関する関数まとめ。

    - base_dir: 出力ベースフォルダ（./output など）
    - progress.xlsx, resume_flag.xlsx のパスを管理
    """

    base_dir: Path

    progress_filename: str = "progress.xlsx"
    resume_flag_filename: str = "resume_flag.xlsx"

    def __post_init__(self) -> None:
        self.base_dir.mkdir(parents=True, exist_ok=True)

    @property
    def progress_path(self) -> Path:
        return self.base_dir / self.progress_filename

    @property
    def resume_flag_path(self) -> Path:
        return self.base_dir / self.resume_flag_filename

    @property
    def result_dir(self) -> Path:
        """
        最終構造をまとめて置く result フォルダ。

        仕様どおり「アプリ直下 ./result」を想定して、
        base_dir=./output の 1 つ上に置く:
            ./output -> ../result
        """
        path = self.base_dir.parent / "result"
        path.mkdir(parents=True, exist_ok=True)
        return path

    # -------- progress.xlsx 関連 --------

    def save_progress(self, rows: List[JobProgressRow]) -> None:
        """全 ID 分の進捗を progress.xlsx に書き出す"""
        df = pd.DataFrame([row.to_dict() for row in rows])
        df.to_excel(self.progress_path, index=False)

    def load_progress(self) -> Optional[List[JobProgressRow]]:
        """
        既存の progress.xlsx があれば読み込み、JobProgressRow のリストとして返す。
        なければ None を返す。
        """
        if not self.progress_path.exists():
            return None

        df = pd.read_excel(self.progress_path)
        rows = [JobProgressRow.from_dict(rec) for _, rec in df.iterrows()]
        return rows

    # -------- resume_flag.xlsx 関連 --------

    def load_resume_flags(self) -> Optional[pd.DataFrame]:
        """resume_flag.xlsx を DataFrame として読み込む（存在しなければ None）"""
        if not self.resume_flag_path.exists():
            return None
        return pd.read_excel(self.resume_flag_path)

    def save_resume_flags(self, df: pd.DataFrame) -> None:
        """resume_flag.xlsx を保存"""
        df.to_excel(self.resume_flag_path, index=False)

    # -------- ヘルパー --------

    def update_progress_row(
        self,
        updated_row: JobProgressRow,
        all_rows: List[JobProgressRow],
    ) -> None:
        """
        all_rows の中から ID をキーにして行を探し、updated_row で置き換えたあと
        progress.xlsx を保存する。
        """
        for i, row in enumerate(all_rows):
            if row.id == updated_row.id:
                all_rows[i] = updated_row
                break
        self.save_progress(all_rows)
