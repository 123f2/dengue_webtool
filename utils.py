"""
utils.py
============================
核心数据处理、查询、可视化函数库
"""

import re
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

GENE_COL = "Gene_Symbol"

DHF_COLOR = "#4C72B0"
DSS_COLOR = "#C44E52"
COMBINED_COLOR = "#55A868"


@st.cache_data(show_spinner="Loading dataset...")
def load_and_preprocess_data(path: str):
    """
    读取原始Excel，完成：
    - GPL201/GPL570 合并重命名为 GSE43777
    - 解析全部患者列名 -> 数据集/严重程度/患者序号
    - 生成数据集元数据汇总表
    """
    df = pd.read_excel(path, sheet_name="Sheet1")

    gpl201_cols = [c for c in df.columns if c.startswith("GPL201_")]
    gpl570_cols = [c for c in df.columns if c.startswith("GPL570_")]

    rename_map = {}
    idx = 1
    for c in gpl201_cols + gpl570_cols:
        rename_map[c] = f"GSE43777_DHFvsDF_{idx}"
        idx += 1
    df = df.rename(columns=rename_map)

    expr_cols = [c for c in df.columns if c != GENE_COL]
    records = []
    for col in expr_cols:
        m = re.match(r"^(GSE\d+)_([A-Za-z]+)vs([A-Za-z]+)_(\d+)$", col)
        if m:
            dataset, severity, control, patient_idx = m.groups()
            records.append({
                "column": col, "dataset": dataset, "severity": severity,
                "control": control, "patient_idx": int(patient_idx),
            })
        else:
            records.append({
                "column": col, "dataset": "UNPARSED", "severity": "UNPARSED",
                "control": "UNPARSED", "patient_idx": None,
            })
    col_meta = pd.DataFrame(records)

    dataset_summary = (
        col_meta.groupby(["dataset", "severity"])
        .size()
        .reset_index(name="n_patients")
        .sort_values(["dataset", "severity"])
        .reset_index(drop=True)
    )

    return df, col_meta, dataset_summary


def query_gene(df: pd.DataFrame, col_meta: pd.DataFrame, gene_symbol: str, mode: str = "combined"):
    """
    查询指定基因在所有患者中的log2FC值。
    mode: "combined" / "DHF" / "DSS"
    """
    match = df[df[GENE_COL].astype(str).str.upper() == gene_symbol.strip().upper()]

    if len(match) == 0:
        return {"found": False, "gene_symbol": gene_symbol, "long_df": None,
                "group_summary": None, "total_valid_n": 0, "total_n": None, "mode": mode}

    row = match.iloc[0]
    real_gene_name = row[GENE_COL]

    expr_cols = col_meta["column"].tolist()
    values = row[expr_cols]

    long_df = col_meta.copy()
    long_df["log2FC"] = long_df["column"].map(values)
    long_df_valid = long_df.dropna(subset=["log2FC"]).copy()

    if mode == "combined":
        filtered = long_df_valid.copy()
        total_n_possible = len(col_meta)
    elif mode in ("DHF", "DSS"):
        filtered = long_df_valid[long_df_valid["severity"] == mode].copy()
        total_n_possible = int((col_meta["severity"] == mode).sum())
    else:
        raise ValueError("mode 必须是 'combined' / 'DHF' / 'DSS'")

    filtered["group"] = filtered["dataset"]
    total_valid_n = len(filtered)

    if total_valid_n == 0:
        group_summary = pd.DataFrame(columns=["group", "n", "mean", "median", "std"])
    else:
        group_summary = (
            filtered.groupby("group")["log2FC"]
            .agg(n="count", mean="mean", median="median", std="std")
            .reset_index()
            .sort_values("group")
            .reset_index(drop=True)
        )

    return {
        "found": True,
        "gene_symbol": real_gene_name,
        "long_df": filtered,
        "group_summary": group_summary,
        "total_valid_n": total_valid_n,
        "total_n": total_n_possible,
        "mode": mode,
    }


def plot_gene_boxplot(query_result: dict):
    if not query_result["found"] or query_result["total_valid_n"] == 0:
        return None

    gene_name = query_result["gene_symbol"]
    long_df = query_result["long_df"].copy()
    total_valid_n = query_result["total_valid_n"]
    total_n = query_result["total_n"]
    mode = query_result["mode"]

    group_n = long_df.groupby("group").size().to_dict()
    long_df["group_label"] = long_df["group"].map(lambda g: f"{g} (n={group_n[g]})")

    order = sorted(long_df["group"].unique())
    order_labels = [f"{g} (n={group_n[g]})" for g in order]

    mode_display = {"combined": "Severe Dengue (DHF+DSS)", "DHF": "DHF", "DSS": "DSS"}[mode]
    color_map = {"combined": COMBINED_COLOR, "DHF": DHF_COLOR, "DSS": DSS_COLOR}

    fig = px.box(
        long_df, x="group_label", y="log2FC", points="all",
        category_orders={"group_label": order_labels},
        title=f"{gene_name} — log2FC Distribution ({mode_display} vs DF)",
        labels={"group_label": "Dataset", "log2FC": "log2 Fold Change"},
        color_discrete_sequence=[color_map[mode]],
    )

    fig.add_hline(y=0, line_dash="dash", line_color="gray", opacity=0.6)

    fig.add_annotation(
        text=f"<b>Overall: n = {total_valid_n} / {total_n} patients analysed</b> "
             f"({total_valid_n/total_n*100:.0f}% detection rate)",
        xref="paper", yref="paper", x=0.5, y=1.10, showarrow=False,
        font=dict(size=13, color="#333333"), align="center",
    )

    fig.update_layout(
        xaxis_tickangle=-30, template="plotly_white", showlegend=False,
        height=560, margin=dict(t=110),
    )
    return fig


def plot_gene_barchart(query_result: dict):
    if not query_result["found"] or query_result["total_valid_n"] == 0:
        return None

    gene_name = query_result["gene_symbol"]
    group_summary = query_result["group_summary"].copy()
    total_valid_n = query_result["total_valid_n"]
    total_n = query_result["total_n"]
    mode = query_result["mode"]

    group_summary["group_label"] = group_summary.apply(
        lambda r: f"{r['group']} (n={int(r['n'])})", axis=1
    )
    group_summary = group_summary.sort_values("group")

    mode_display = {"combined": "Severe Dengue (DHF+DSS)", "DHF": "DHF", "DSS": "DSS"}[mode]
    color_map = {"combined": COMBINED_COLOR, "DHF": DHF_COLOR, "DSS": DSS_COLOR}

    fig = px.bar(
        group_summary, x="group_label", y="mean", error_y="std",
        title=f"{gene_name} — Mean log2FC by Dataset ({mode_display} vs DF)",
        labels={"group_label": "Dataset", "mean": "Mean log2 Fold Change"},
        color_discrete_sequence=[color_map[mode]],
    )

    fig.add_hline(y=0, line_dash="dash", line_color="gray", opacity=0.6)

    fig.add_annotation(
        text=f"<b>Overall: n = {total_valid_n} / {total_n} patients analysed</b> "
             f"({total_valid_n/total_n*100:.0f}% detection rate)",
        xref="paper", yref="paper", x=0.5, y=1.10, showarrow=False,
        font=dict(size=13, color="#333333"), align="center",
    )

    fig.update_layout(
        xaxis_tickangle=-30, template="plotly_white", showlegend=False,
        height=520, margin=dict(t=110),
    )
    return fig

# ============ 追加到 utils.py 末尾 ============

@st.cache_data(show_spinner="Computing overall distribution across datasets...")
def get_overall_dataset_distribution(df: pd.DataFrame, col_meta: pd.DataFrame, sample_per_dataset: int = 20000):
    expr_cols = col_meta["column"].tolist()
    melted = df[expr_cols].melt(var_name="column", value_name="log2FC").dropna()
    melted = melted.merge(col_meta[["column", "dataset", "severity"]], on="column", how="left")

    sampled = melted.groupby("dataset", group_keys=False).apply(
        lambda x: x.sample(n=min(len(x), sample_per_dataset), random_state=42)
    )
    return sampled


def plot_overall_distribution(melted_df: pd.DataFrame):
    """
    整体分布箱线图：x轴=7个数据集，y轴=该数据集下所有基因所有患者的log2FC值。
    注意：这里用 points=False（不画散点），因为数据量是百万级别，画散点会导致
    浏览器卡顿；这张图的目的是看整体分布形状，不是看单个数据点。
    """
    order = sorted(melted_df["dataset"].unique())
    n_values_per_dataset = melted_df.groupby("dataset").size().to_dict()

    melted_df = melted_df.copy()
    melted_df["dataset_label"] = melted_df["dataset"].map(
        lambda d: f"{d} ({n_values_per_dataset[d]:,} values)"
    )
    order_labels = [f"{d} ({n_values_per_dataset[d]:,} values)" for d in order]

    fig = px.box(
        melted_df,
        x="dataset_label",
        y="log2FC",
        points=False,
        category_orders={"dataset_label": order_labels},
        title="Overall log2FC Distribution per Dataset (all genes combined)",
        labels={"dataset_label": "Dataset", "log2FC": "log2 Fold Change"},
        color_discrete_sequence=["#4C72B0"],
    )
    fig.add_hline(y=0, line_dash="dash", line_color="gray", opacity=0.6)
    fig.update_layout(
        xaxis_tickangle=-30,
        template="plotly_white",
        showlegend=False,
        height=520,
    )
    return fig
    
