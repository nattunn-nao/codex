from __future__ import annotations

import streamlit as st

from conformer_app.core.config import (
    EmbedConfig,
    MMConfig,
    UMAConfig,
    QMConfig,
    OutputConfig,
)


# =========================
# mol_gen（コンフォマー生成） UI
# =========================

def ui_embed_block(default: EmbedConfig) -> EmbedConfig:
    """
    コンフォマー生成（ETKDG）と、塩・イオンペアなど
    '.' で分割される複数フラグメントの空間配置パラメータを設定する。
    """
    with st.expander("コンフォマー生成", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            num_confs = st.number_input(
                "最大コンフォマー数",
                min_value=1,
                max_value=5000,
                value=int(default.num_confs),
                step=100,
            )
            prune = st.number_input(
                "生成時重複除去RMS (Å, 負値ならオフ)",
                min_value=-1.0,
                max_value=5.0,
                value=float(default.prune_rms),
                step=0.1,
            )
        with c2:
            etkdg_version = st.selectbox(
                "ETKDGバージョン",
                options=["ETKDGv3", "ETKDGv2", "ETKDG"],
                index={"ETKDGv3": 0, "ETKDGv2": 1, "ETKDG": 2}.get(default.etkdg_version, 0),
            )
            seed = st.number_input(
                "乱数シード（-1: ランダム）",
                min_value=-1,
                max_value=10_000_000,
                value=int(default.random_seed),
                step=1,
            )
        with c3:
            small_ring = st.checkbox(
                "Small ring torsions",
                value=bool(default.small_ring),
            )
            macrocycle = st.checkbox(
                "Macrocycle torsions",
                value=bool(default.macrocycle),
            )

        st.markdown("**複数フラグメント（塩・イオンペア）の配置**")
        g1, g2 = st.columns(2)
        with g1:
            frag_gap_min = st.number_input(
                "フラグメント間の最小表面距離 [Å]",
                min_value=0.0,
                max_value=20.0,
                value=float(getattr(default, "frag_gap_min", 2.0)),
                step=0.1,
            )
        with g2:
            frag_gap_max = st.number_input(
                "フラグメント間の最大表面距離 [Å]",
                min_value=0.0,
                max_value=50.0,
                value=float(getattr(default, "frag_gap_max", 6.0)),
                step=0.1,
            )

        st.caption(
            "SMILES に [.] を含む場合\n"
            "最も重いフラグメントを中心に置き、そのvdw表面から他フラグメントの“表面”までの\n"
            "距離が [最小, 最大] の範囲に入るように配置します。"
        )

    return EmbedConfig(
        num_confs=int(num_confs),
        prune_rms=float(prune),
        etkdg_version=str(etkdg_version),
        random_seed=int(seed),
        small_ring=bool(small_ring),
        macrocycle=bool(macrocycle),
        frag_gap_min=float(frag_gap_min),
        frag_gap_max=float(frag_gap_max),
    )


# =========================
# MM 設定 UI
# =========================

def ui_mm_block(default: MMConfig) -> MMConfig:
    """
    力場 / 反復数に加え、エネルギーによる重複除去と
    UMA に渡す構造数（エネルギー上位 N）の設定。
    """
    with st.expander("MM 最適化 / 構造選抜", expanded=False):
        c1, c2 = st.columns(2)
        with c1:
            ff_mode = st.selectbox(
                "力場",
                options=["MMFF94s->UFF", "MMFF94", "UFF"],
                index={"MMFF94s->UFF": 0, "MMFF94": 1, "UFF": 2}.get(default.ff_mode, 0),
            )
        with c2:
            max_iters = st.number_input(
                "反復回数（MM 最適化）",
                min_value=50,
                max_value=5000,
                value=int(default.max_iters),
                step=50,
            )

        st.markdown("**エネルギーによる重複除去 & 窓**")
        c3, c4 = st.columns(2)
        with c3:
            equal_tol = st.number_input(
                "同値判定 ΔE (kcal/mol)",
                min_value=0.0,
                max_value=0.1,
                value=float(getattr(default, "equal_tol", 0.003)),
                step=0.001,
                format="%.4f",
            )
        with c4:
            energy_window = st.number_input(
                "エネルギー窓 (E ≤ Emin + 窓, kcal/mol)",
                min_value=0.0,
                max_value=100.0,
                value=float(getattr(default, "energy_window", 10.0)),
                step=1.0,
            )

        st.markdown("**MMopt 後に UMA に渡す構造数**")
        c5, c6 = st.columns(2)
        with c5:
            max_keep = st.number_input(
                "UMA に渡す最大構造数",
                min_value=1,
                max_value=200,
                value=int(getattr(default, "max_keep", 10)),
                step=1,
            )

    return MMConfig(
        ff_mode=str(ff_mode),
        max_iters=int(max_iters),
        equal_tol=float(equal_tol),
        energy_window=float(energy_window),
        rmsd_min=getattr(default, "rmsd_min", 0.8),
        max_keep=int(max_keep),
    )


# =========================
# UMA 設定 UI
# =========================

def ui_uma_block(default: UMAConfig) -> tuple[UMAConfig, bool]:
    """
    UMA ポテンシャル + LBFGS の設定。

    戻り値:
        (uma_config, use_uma_flag)

        use_uma_flag が False の場合、
        app.py 側で UMA 最適化をスキップし、
        post_MM から直接 QM 構造を選ぶモードに切り替える。
    """
    with st.expander("UMA 設定", expanded=False):
        # UMA を使う / 使わないのトグル
        use_uma = st.checkbox(
            "UMA ポテンシャルで最適化を行う",
            value=True,
            help="チェックを外すと UMA 計算をスキップし、MM の結果から直接 QM 構造を選びます。",
        )

        model_name = st.text_input("モデル", value=default.model_name)

        c1, c2 = st.columns(2)
        with c1:
            fmax = st.number_input(
                "LBFGS 収束閾値 (eV/Å)",
                min_value=0.0,
                max_value=1.0,
                value=float(default.fmax),
                step=0.01,
            )
        with c2:
            max_steps = st.number_input(
                "LBFGS 最大反復",
                min_value=10,
                max_value=100_000,
                value=int(default.max_steps),
                step=10,
            )

        c3, c4 = st.columns(2)
        with c3:
            device = st.selectbox(
                "デバイス",
                options=["cpu", "cuda"],
                index=0 if default.device == "cpu" else 1,
            )
        with c4:
            debug = st.checkbox(
                "デバッグ出力",
                value=bool(default.debug),
            )

        st.markdown("系の電荷 / 多重度")
        cc, cm = st.columns(2)
        with cc:
            charge = st.number_input("電荷 (UMA charge)", -10, 10, int(default.charge), 1)
        with cm:
            spin = st.number_input("スピン多重度 (2S+1)", 1, 10, int(default.spin), 1)

    uma_cfg = UMAConfig(
        model_name=str(model_name),
        device=str(device),
        fmax=float(fmax),
        max_steps=int(max_steps),
        charge=int(charge),
        spin=int(spin),
        task_name=default.task_name,
        debug=bool(debug),
    )
    return uma_cfg, bool(use_uma)


# =========================
# GAMESS / QM 設定 UI
# =========================

def ui_gamess_block(default: QMConfig) -> QMConfig:
    """
    GAMESS のテンプレート / 実行スクリプトおよび CPU リソースを設定。
    ICHARG は通常、Progress.xlsx に書かれた Formal_Charge が優先され、
    ここでの値は「うまく取れなかったときのフォールバック」として使われる。
    """
    with st.expander("GAMESS / QM 設定", expanded=False):
        st.markdown("**テンプレート / スクリプト**")
        tpl = st.text_input("gms.inp パス", value=default.template_path)
        bat = st.text_input("run.bat パス", value=default.run_bat_path)
        ps1 = st.text_input("gms_progress.ps1 パス", value=default.ps1_path)

        st.markdown("**計算資源**")
        c1, c2 = st.columns(2)
        with c1:
            ncpus = st.number_input(
                "NCPUS（ジョブあたり）",
                min_value=1,
                max_value=128,
                value=int(default.ncpus),
                step=1,
            )
        with c2:
            maxp = st.number_input(
                "同時実行ジョブ数",
                min_value=1,
                max_value=32,
                value=int(default.max_parallel),
                step=1,
            )

        st.markdown("**系の電荷 / 多重度**")
        c3, c4 = st.columns(2)
        with c3:
            icharge = st.number_input(
                "電荷 (ICHARG, fallback)",
                min_value=-10,
                max_value=10,
                value=int(default.icharge),
                step=1,
            )
        with c4:
            mult = st.number_input(
                "多重度 (MULT)",
                min_value=1,
                max_value=10,
                value=int(default.mult),
                step=1,
            )

        pointgroup = st.text_input("ポイントグループ (%POINTGROUP%)", value=default.pointgroup)
        smiles_comment = st.text_input("コメント（%COMMENT_SMILES_N%）", value=default.smiles_comment)

    return QMConfig(
        template_path=str(tpl).strip(),
        run_bat_path=str(bat).strip(),
        ps1_path=str(ps1).strip(),
        ncpus=int(ncpus),
        max_parallel=int(maxp),
        icharge=int(icharge),
        mult=int(mult),
        pointgroup=str(pointgroup),
        smiles_comment=str(smiles_comment),
    )


# =========================
# 出力設定 UI
# =========================

def ui_output_block(default: OutputConfig) -> OutputConfig:
    with st.expander("出力設定", expanded=False):
        out_base = st.text_input("出力ベースフォルダ", value=default.base_dir)
        move_com = st.checkbox(
            "XYZ 出力時に重心を原点へ移動",
            value=bool(default.move_to_com_origin),
        )
    return OutputConfig(
        base_dir=str(out_base),
        move_to_com_origin=bool(move_com),
    )
