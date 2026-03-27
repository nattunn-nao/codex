from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from .config import load_default_config
from .registry import load_surface_registry
from .session_state import init_session_state, push_log
from .surface_loader import load_surface
from .workflow import run_batch


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    st.set_page_config(page_title="CMP吸着構造探索", layout="wide")
    st.title("CMP吸着構造探索アプリ")

    init_session_state()

    cfg = load_default_config(root)
    entries = load_surface_registry(root)
    labels = {f"{e.label} ({e.sid})": e for e in entries}

    st.sidebar.header("計算条件")
    selected_label = st.sidebar.selectbox("表面モデル", options=list(labels.keys()))
    selected = labels[selected_label]

    pose = cfg.data["pose"]
    pose["xy_step"] = st.sidebar.number_input("xyステップ [Å]", min_value=0.5, value=float(pose["xy_step"]), step=0.1)
    pose["max_poses"] = st.sidebar.number_input("最大初期構造数", min_value=1, value=int(pose["max_poses"]), step=1)

    z_input = st.sidebar.text_input("高さリスト z [Å]", value=",".join(map(str, pose["z_values"])))
    r_input = st.sidebar.text_input("回転角 [deg]", value=",".join(map(str, pose["rotations_deg"])))
    pose["z_values"] = [float(x) for x in z_input.split(",") if x.strip()]
    pose["rotations_deg"] = [float(x) for x in r_input.split(",") if x.strip()]

    st.write(f"入力分子フォルダ: `{cfg.input_adsorbates_dir}` 内の `.xyz` を逐次処理します。")
    if st.button("探索を開始"):
        surface = load_surface(selected.sid, selected.path)
        push_log(f"Surface loaded: {selected.sid}, fixed={len(surface.fixed_indices)}")
        summaries = run_batch(surface=surface, cfg=cfg, root=root)
        rows = [s.__dict__ for s in summaries]
        df = pd.DataFrame(rows)
        st.subheader("実行サマリー")
        st.dataframe(df, use_container_width=True)
        for _, row in df.iterrows():
            push_log(f"{row['adsorbate_file']}: best={row['best_energy']:.4f} eV, report={row['report_csv']}")

    st.subheader("ログ")
    st.code("\n".join(st.session_state["logs"]) or "(no logs)")


if __name__ == "__main__":
    main()
