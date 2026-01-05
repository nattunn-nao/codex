from __future__ import annotations

"""
mm_layer.builder

複数フラグメントを含む SMILES（例: "O=S(=O)([O-])[O-].[NH4+].[NH4+]"} に対して、

- フラグメントごとに ETKDG で多数のコンフォマーを生成
- 各フラグメントの重心を原点に移動し「半径」を定義
- 一番重いフラグメントを原点、他フラグメントをシェル上にランダム配置
- そのクラスター構造を RDKit Mol の複数 Conformer として構築

を行う。

公開関数:
    build_clustered_multiconformer_mol(smiles, embed_cfg) -> (mol, conf_ids)
"""

import math
import random
from typing import List, Dict, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

from conformer_app.core.config import EmbedConfig


# ==============================
# 距離・シェル配置に関するデフォルト
# ==============================

# 代表フラグメント表面から見たサブフラグメント最近接端までの距離 [Å]
DEFAULT_SURFACE_GAP_MIN = 1.5   # 近すぎない下限
DEFAULT_SURFACE_GAP_MAX = 3.0   # 離れすぎない上限


# ==============================
# ユーティリティ
# ==============================

def _split_smiles(smiles: str) -> List[str]:
    """'.' で SMILES をフラグメントに分割する。"""
    parts = [p for p in smiles.replace(" ", "").split(".") if p]
    if not parts:
        raise ValueError(f"フラグメントが見つかりません: {smiles}")
    return parts


def _count_heavy_atoms(smiles: str) -> int:
    """重原子数（C, N, O, ...）を数える。"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES 読み込み失敗: {smiles}")
    return mol.GetNumHeavyAtoms()


def _compute_radius(coords: List[Tuple[str, float, float, float]]) -> float:
    """
    coords: [(sym, x, y, z), ...]（重心が原点の前提）
    原点から最も遠い点までの距離を「半径」とみなす。
    """
    max_r2 = 0.0
    for _, x, y, z in coords:
        r2 = x * x + y * y + z * z
        if r2 > max_r2:
            max_r2 = r2
    return math.sqrt(max_r2)


def _random_unit_vector() -> Tuple[float, float, float]:
    """3次元の一様乱数方向ベクトル。"""
    z = random.uniform(-1.0, 1.0)
    t = random.uniform(0.0, 2.0 * math.pi)
    r_xy = math.sqrt(max(0.0, 1.0 - z * z))
    x = r_xy * math.cos(t)
    y = r_xy * math.sin(t)
    return x, y, z


def _dist(p: Tuple[float, float, float], q: Tuple[float, float, float]) -> float:
    """2点間距離。"""
    dx = p[0] - q[0]
    dy = p[1] - q[1]
    dz = p[2] - q[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _get_etkdg_params(embed_cfg: EmbedConfig):
    """
    EmbedConfig に基づいて ETKDG パラメータを構築。
    """
    ver = (embed_cfg.etkdg_version or "ETKDGv3").upper()
    if ver == "ETKDGV2":
        params = AllChem.ETKDGv2()
    elif ver == "ETKDG":
        params = AllChem.ETKDG()
    else:
        params = AllChem.ETKDGv3()

    # 生成時重複除去 RMS
    if embed_cfg.prune_rms is not None and embed_cfg.prune_rms >= 0.0:
        params.pruneRmsThresh = float(embed_cfg.prune_rms)
    else:
        params.pruneRmsThresh = -1.0  # 無効化

    # 乱数シード
    if embed_cfg.random_seed is not None and embed_cfg.random_seed >= 0:
        params.randomSeed = int(embed_cfg.random_seed)
    else:
        params.randomSeed = -1  # 完全ランダム

    # small ring / macrocycle
    try:
        params.useSmallRingTorsions = bool(embed_cfg.small_ring)
    except Exception:
        pass
    try:
        params.useMacrocycleTorsions = bool(embed_cfg.macrocycle)
    except Exception:
        pass

    return params


# ==============================
# コンフォマー生成（フラグメントごと）
# ==============================

def _build_fragment_conformers_from_mol(
    mol: Chem.Mol,
    n_confs: int,
    embed_cfg: EmbedConfig,
) -> List[Dict[str, object]]:
    """
    1つの「H 追加済みフラグメント mol」に対して、
    最大 n_confs 個のコンフォマーを生成し、
    重心を原点に平行移動した座標と「半径」を返す。

    戻り値: list of
        {
          "coords": [(sym, x, y, z), ...],  # 重心が原点
          "radius": float
        }
    """
    # clone しておく（外側の mol に影響させたくない場合）
    mol = Chem.Mol(mol)

    params = _get_etkdg_params(embed_cfg)

    conf_ids = AllChem.EmbedMultipleConfs(
        mol,
        numConfs=int(n_confs),
        params=params,
    )
    if not conf_ids:
        raise RuntimeError("コンフォマー埋め込みに失敗しました。")

    # ざっくり UFF で緩めに最適化（テンプレの形を整える程度）
    AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=0)

    conformers: List[Dict[str, object]] = []
    n_atoms = mol.GetNumAtoms()

    for cid in conf_ids:
        conf = mol.GetConformer(cid)

        raw: List[Tuple[str, float, float, float]] = []
        for i in range(n_atoms):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            raw.append((atom.GetSymbol(), float(pos.x), float(pos.y), float(pos.z)))

        # 重心を原点に移動
        cx = sum(x for _, x, _, _ in raw) / n_atoms
        cy = sum(y for _, _, y, _ in raw) / n_atoms
        cz = sum(z for _, _, _, z in raw) / n_atoms

        centered: List[Tuple[str, float, float, float]] = []
        for sym, x, y, z in raw:
            centered.append((sym, x - cx, y - cy, z - cz))

        radius = _compute_radius(centered)
        conformers.append({"coords": centered, "radius": radius})

    return conformers


# ==============================
# shell 内にサブ構造を配置
# ==============================

def _place_fragments_on_shell(
    ordered_frags: List[Dict[str, object]],
    gap_min: float,
    gap_max: float,
) -> List[Tuple[float, float, float]]:
    """
    ordered_frags: [frag0, frag1, ...]
      frag = {"radius": r, ...}

    0番目のフラグメントを原点に置き、
    代表構造からの距離に基づいて shell を定義：

      r_main = radius(main)
      R_inner = r_main + gap_min
      R_outer = r_main + gap_max

    各サブ構造 i について：
      - サブ構造の最近接末端の半径 r_surface を
        [R_inner, R_outer] から乱数で選ぶ
      - 中心距離 d = r_surface + r_i として、ランダムな方向に配置
      - 他フラグメントと重ならないように
        d_ij >= r_i + r_j + gap_min をチェック
    """
    if gap_max <= gap_min:
        raise ValueError("gap_max は gap_min より大きくしてください。")

    n = len(ordered_frags)
    if n == 0:
        return []

    radii = [float(f["radius"]) for f in ordered_frags]
    r_main = radii[0]

    R_inner = r_main + gap_min
    R_outer = r_main + gap_max

    placements: List[Tuple[float, float, float]] = [(0.0, 0.0, 0.0)]  # 代表は原点

    if n == 1:
        return placements

    max_trials = 300

    for i in range(1, n):
        r_i = radii[i]
        placed = False

        for _ in range(max_trials):
            # shell 内の「最近接末端半径」をランダムに選ぶ
            r_surface = random.uniform(R_inner, R_outer)
            # 中心距離 = 末端半径 + 自身の半径
            d = r_surface + r_i

            ux, uy, uz = _random_unit_vector()
            cand = (d * ux, d * uy, d * uz)

            ok = True
            # 他フラグメントとの重なりをチェック（最小距離のみ）
            for j in range(i):
                r_j = radii[j]
                d_ij = _dist(cand, placements[j])
                d_min = r_i + r_j + gap_min
                if d_ij < d_min:
                    ok = False
                    break

            if ok:
                placements.append(cand)
                placed = True
                break

        if not placed:
            # 条件が厳しすぎたら、外側に少し逃がす
            ux, uy, uz = _random_unit_vector()
            d = R_outer + r_i + gap_min
            placements.append((d * ux, d * uy, d * uz))

    return placements


# ==============================
# メイン: RDKit Mol + Conformers 構築
# ==============================

def build_clustered_multiconformer_mol(
    smiles: str,
    embed_cfg: EmbedConfig,
) -> Tuple[Chem.Mol, List[int]]:
    """
    1つの（複数フラグメントを含みうる）SMILES に対して:

    - フラグメントごとに多数のコンフォマーを生成
    - 一番重いフラグメントを代表構造にし、他フラグメントをシェル上に配置
    - そのクラスターを RDKit Mol の各 Conformer として構築

    Parameters
    ----------
    smiles : str
        入力 SMILES（例: "O=S(=O)([O-])[O-].[NH4+].[NH4+]"）
    embed_cfg : EmbedConfig
        コンフォマー生成に関する設定（num_confs 等）。

    Returns
    -------
    mol : Chem.Mol
        複数フラグメントを含む RDKit 分子（H 付き）。
    conf_ids : list of int
        生成された Conformer の ID 一覧。
    """
    # Python の乱数もシード可能にしておく（再現性がほしい場合用）
    if embed_cfg.random_seed is not None and embed_cfg.random_seed >= 0:
        random.seed(int(embed_cfg.random_seed))

    frag_smis = _split_smiles(smiles)

    # コンフォマープールのサイズ
    # →「生成したいクラスター構造数」と同程度か、それ以上
    num_structures = int(embed_cfg.num_confs)
    conf_pool_size = max(num_structures, 50)

    # ---- フラグメントごとにテンプレ構造（コンフォマー群）を準備 ----
    frag_templates: List[Dict[str, object]] = []
    frag_mols: List[Chem.Mol] = []

    for fs in frag_smis:
        heavy = _count_heavy_atoms(fs)

        base = Chem.MolFromSmiles(fs)
        if base is None:
            raise ValueError(f"フラグメント SMILES 読み込みに失敗: {fs}")

        molH = Chem.AddHs(base)
        # H 追加後のフラグメントを明示的にサニタイズして RingInfo 等を初期化
        Chem.SanitizeMol(molH)
        molH.UpdatePropertyCache(strict=False)

        frag_mols.append(molH)

        confs = _build_fragment_conformers_from_mol(
            mol=molH,
            n_confs=conf_pool_size,
            embed_cfg=embed_cfg,
        )
        if not confs:
            raise RuntimeError(f"フラグメント {fs} のコンフォマー生成に失敗")

        frag_templates.append(
            {
                "smiles": fs,
                "heavy": heavy,
                "confs": confs,
            }
        )

    # ---- フラグメントを 1 つの Mol に結合 ----
    if not frag_mols:
        raise RuntimeError("フラグメント分子が 1 つも生成できませんでした。")

    combined = frag_mols[0]
    for m in frag_mols[1:]:
        combined = Chem.CombineMols(combined, m)

    rw = Chem.RWMol(combined)
    Chem.SanitizeMol(rw)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)

    # フラグメントごとの atom index リストを作る
    frag_atom_lists: List[List[int]] = []
    offset = 0
    for m in frag_mols:
        n = m.GetNumAtoms()
        frag_atom_lists.append(list(range(offset, offset + n)))
        offset += n

    if offset != mol.GetNumAtoms():
        raise RuntimeError(
            f"フラグメントの原子数 ({offset}) と結合後 Mol の原子数 ({mol.GetNumAtoms()}) が一致しません。"
        )

    # ---- クラスター構造を num_structures 個作成し、Conformer として追加 ----
    gap_min = DEFAULT_SURFACE_GAP_MIN
    gap_max = DEFAULT_SURFACE_GAP_MAX

    conf_ids: List[int] = []

    # 代表構造（main fragment）は「重原子数が最大」のもの
    main_frag_index = max(range(len(frag_templates)), key=lambda i: frag_templates[i]["heavy"])

    for _ in range(num_structures):
        # 各フラグメントからコンフォマーをランダムに 1 つ選ぶ
        chosen_confs: List[Dict[str, object]] = []
        for tmpl in frag_templates:
            conf = random.choice(tmpl["confs"])  # {"coords": [...], "radius": ...}
            chosen_confs.append(conf)

        # 代表構造を先頭にして並べ替え
        ordered_indices = [main_frag_index] + [
            i for i in range(len(frag_templates)) if i != main_frag_index
        ]
        ordered_confs = [chosen_confs[i] for i in ordered_indices]

        ordered_for_shell = [{"radius": c["radius"]} for c in ordered_confs]
        offsets = _place_fragments_on_shell(
            ordered_for_shell,
            gap_min=gap_min,
            gap_max=gap_max,
        )

        # RDKit Conformer を構築
        conf = Chem.Conformer(mol.GetNumAtoms())

        # ordered_indices / offsets を元に、元のフラグメント index にマッピング
        for ord_idx, frag_idx in enumerate(ordered_indices):
            dx, dy, dz = offsets[ord_idx]
            coords = chosen_confs[frag_idx]["coords"]  # [(sym, x, y, z), ...]
            atom_indices = frag_atom_lists[frag_idx]

            if len(coords) != len(atom_indices):
                raise RuntimeError(
                    f"フラグメント {frag_idx} の座標数 ({len(coords)}) と "
                    f"原子数 ({len(atom_indices)}) が一致しません。"
                )

            for atom_idx, (_, x, y, z) in zip(atom_indices, coords):
                conf.SetAtomPosition(int(atom_idx), (float(x + dx), float(y + dy), float(z + dz)))

        cid = mol.AddConformer(conf, assignId=True)
        conf_ids.append(cid)

    return mol, conf_ids
