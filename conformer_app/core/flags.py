from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ControlFlags:
    """
    実行中の挙動を制御するフラグ。

    - stop_immediately: 現在のステップ終了後に即座に全体を止める
    - stop_after_current: 現在の ID について QM まで終わったタイミングで止める
    """

    stop_immediately: bool = False
    stop_after_current: bool = False

    @property
    def should_stop_immediately(self) -> bool:
        return self.stop_immediately

    @property
    def should_stop_after_current(self) -> bool:
        return self.stop_after_current

    def clear(self) -> None:
        """フラグをリセット"""
        self.stop_immediately = False
        self.stop_after_current = False
