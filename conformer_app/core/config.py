from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class OutputConfig:
    """
    出力関連の設定
    """
    base_dir: str = "./output"          # ID ごとのフォルダを掘るルート
    move_to_com_origin: bool = True     # XYZ 出力時に重心を原点に移動するか


@dataclass
class EmbedConfig:
    # ETKDG / コンフォマー生成
    num_confs: int = 1000
    prune_rms: float = 0.3
    etkdg_version: str = "ETKDGv3"
    random_seed: int = 251027
    small_ring: bool = True
    macrocycle: bool = True
    frag_gap_min: float = 1.0
    frag_gap_max: float = 2.0


@dataclass
class MMConfig:
    """
    MM 最適化関連
    """
    ff_mode: str = "MMFF94s->UFF"
    max_iters: int = 500
    equal_tol: float = 0.003
    energy_window: float = 10.0
    rmsd_min: float = 0.8
    max_keep: int = 10   # ← UMA に渡す最大構造数



@dataclass
class UMAConfig:
    """
    UMA ポテンシャル関連の設定
    """
    model_name: str = "uma-s-1p1"       # fairchem のチェックポイント名
    device: str = "cpu"                 # "cpu" or "cuda"
    fmax: float = 0.01                  # 収束閾値 [eV/Å]
    max_steps: int = 500                # LBFGS 最大ステップ数

    # 将来のプロトン化状態対応用
    charge: int = 0
    spin: int = 1

    task_name: str = "omol"
    debug: bool = False


@dataclass
class QMConfig:
    """
    GAMESS/QM 関連の設定
    """
    # デフォルトで ./template 配下を参照する
    template_path: str = "./template/gms.inp"          # gms.inp
    run_bat_path: str = "./template/run.bat"           # run.bat
    ps1_path: str = "./template/gms_progress.ps1"      # gms_progress.ps1

    # 計算資源
    ncpus: int = 8                      # 1 ジョブあたりの CPU コア数
    max_parallel: int = 3               # 同時実行ジョブ数

    # 系の電荷 / 多重度
    icharge: int = 0                    # ICHARG
    mult: int = 1                       # MULT

    # その他
    pointgroup: str = "C1"
    smiles_comment: str = ""            # %COMMENT_SMILES_N% に入れる文字列


@dataclass
class AppConfig:
    """
    アプリ全体の設定ラッパー
    """
    output: OutputConfig = field(default_factory=OutputConfig)
    embed:  EmbedConfig  = field(default_factory=EmbedConfig)
    mm:     MMConfig     = field(default_factory=MMConfig)
    uma:    UMAConfig    = field(default_factory=UMAConfig)
    qm:     QMConfig     = field(default_factory=QMConfig)
