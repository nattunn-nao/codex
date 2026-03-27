from __future__ import annotations

import streamlit as st


def init_session_state() -> None:
    st.session_state.setdefault("logs", [])


def push_log(message: str) -> None:
    st.session_state["logs"].append(message)
