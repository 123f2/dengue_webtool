"""
utils.py
============================
核心数据处理、查询、可视化函数库
"""

import re
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import pearsonr
import streamlit as st

GENE_COL = "Gene_Symbol"

DHF_COLOR = "#4C72B0"
DSS_COLOR = "#C44E52"
COMBINED_COLOR = "#55A868"
UP_COLOR = "#C44E52"
DOWN_COLOR = "#4C72B0"


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

# =========================================================
# 4.1 数据集内 Top 基因排名 (全景散点图 + 表格)
# =========================================================
@st.cache_data(show_spinner="Computing gene landscape for this dataset...")
def get_gene_landscape(df: pd.DataFrame, col_meta: pd.DataFrame, dataset_name: str,
                        mode: str = "combined", min_patients: int = 3):
    """
    返回该数据集下全部基因的地形图数据，按均值log2FC排序。
    x_position 用基因在原始df里的行号(跟数值无关，避免高亮点在绘图时挤在一起)
    """
    if mode == "combined":
        cols_meta_f = col_meta[col_meta["dataset"] == dataset_name]
    else:
        cols_meta_f = col_meta[(col_meta["dataset"] == dataset_name) & (col_meta["severity"] == mode)]

    cols = cols_meta_f["column"].tolist()
    if len(cols) == 0:
        return None

    sub = df[cols]
    means = sub.mean(axis=1, skipna=True)
    counts = sub.notna().sum(axis=1)

    result = pd.DataFrame({
        "gene": df[GENE_COL],
        "mean_log2FC": means,
        "n_patients": counts,
        "x_position": df.index,
    })
    result = result[result["n_patients"] >= min_patients].copy()
    result = result.sort_values("mean_log2FC", ascending=True).reset_index(drop=True)
    result["rank"] = result.index + 1

    return result


def mark_top_hits(landscape_df: pd.DataFrame, top_n: int = 10):
    """在地形图数据两端各标记top_n个基因为高亮 (Up / Down)"""
    result = landscape_df.copy()
    result["highlight"] = "No"
    result.loc[result.index[:top_n], "highlight"] = "Down"
    result.loc[result.index[-top_n:], "highlight"] = "Up"
    return result


def plot_gene_rank_scatter(marked_df: pd.DataFrame, dataset_name: str, mode_label: str = "combined"):
    """GeneRaMeN风格全景散点图：灰色背景全部基因 + 两端高亮上色(不标常驻文字，靠hover)"""
    normal_df = marked_df[marked_df["highlight"] == "No"]
    up_df = marked_df[marked_df["highlight"] == "Up"]
    down_df = marked_df[marked_df["highlight"] == "Down"]

    fig = go.Figure()

    fig.add_trace(go.Scattergl(
        x=normal_df["x_position"], y=normal_df["mean_log2FC"],
        mode="markers",
        marker=dict(color="lightgray", size=5, opacity=0.5),
        name="Other genes",
        hoverinfo="skip",
    ))

    fig.add_trace(go.Scatter(
        x=down_df["x_position"], y=down_df["mean_log2FC"],
        mode="markers",
        marker=dict(color=DOWN_COLOR, size=9),
        name="Top Downregulated",
        customdata=down_df["gene"],
        hovertemplate="<b>%{customdata}</b><br>log2FC=%{y:.3f}<extra></extra>",
    ))

    fig.add_trace(go.Scatter(
        x=up_df["x_position"], y=up_df["mean_log2FC"],
        mode="markers",
        marker=dict(color=UP_COLOR, size=9),
        name="Top Upregulated",
        customdata=up_df["gene"],
        hovertemplate="<b>%{customdata}</b><br>log2FC=%{y:.3f}<extra></extra>",
    ))

    fig.add_hline(y=0, line_dash="dash", line_color="gray", opacity=0.5)

    fig.update_layout(
        title=f"{dataset_name} — Gene Ranking Landscape ({mode_label})",
        xaxis_title="Gene index (original order, not sorted by value)",
        yaxis_title="Mean log2FC",
        template="plotly_white",
        height=550,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return fig


def get_top_hits_table(marked_df: pd.DataFrame):
    """从mark_top_hits()结果提取Up/Down两张Top基因表格"""
    hits = marked_df[marked_df["highlight"] != "No"].copy()

    up_table = hits[hits["highlight"] == "Up"].sort_values("mean_log2FC", ascending=False).reset_index(drop=True)
    down_table = hits[hits["highlight"] == "Down"].sort_values("mean_log2FC", ascending=True).reset_index(drop=True)

    up_table.insert(0, "Rank", range(1, len(up_table) + 1))
    down_table.insert(0, "Rank", range(1, len(down_table) + 1))

    up_table = up_table[["Rank", "gene", "mean_log2FC", "n_patients"]].rename(
        columns={"gene": "Gene", "mean_log2FC": "Mean log2FC", "n_patients": "N Patients"}
    )
    down_table = down_table[["Rank", "gene", "mean_log2FC", "n_patients"]].rename(
        columns={"gene": "Gene", "mean_log2FC": "Mean log2FC", "n_patients": "N Patients"}
    )

    return up_table, down_table


# =========================================================
# 4.2 基因跨数据集排名查询
# =========================================================
def gene_rank_across_datasets(df: pd.DataFrame, col_meta: pd.DataFrame, gene_symbol: str,
                                mode: str = "combined", min_patients: int = 3):
    """对指定基因，计算它在每个数据集里的：均值log2FC、排名、百分位数、检出患者数"""
    match = df[df[GENE_COL].astype(str).str.upper() == gene_symbol.strip().upper()]
    if len(match) == 0:
        return None, gene_symbol

    real_name = match.iloc[0][GENE_COL]
    gene_row_idx = match.index[0]

    records = []
    all_datasets = sorted(col_meta["dataset"].unique())
    for dataset_name in all_datasets:
        if mode == "combined":
            cols_meta_f = col_meta[col_meta["dataset"] == dataset_name]
        else:
            cols_meta_f = col_meta[(col_meta["dataset"] == dataset_name) & (col_meta["severity"] == mode)]

        cols = cols_meta_f["column"].tolist()
        if len(cols) == 0:
            continue

        sub = df[cols]
        means = sub.mean(axis=1, skipna=True)
        counts = sub.notna().sum(axis=1)

        gene_n_patients = int(counts.loc[gene_row_idx])
        gene_mean = means.loc[gene_row_idx]

        if gene_n_patients < min_patients or pd.isna(gene_mean):
            records.append({
                "dataset": dataset_name, "mean_log2FC": np.nan,
                "rank": None, "n_valid_genes": None, "percentile": None,
                "n_patients": gene_n_patients,
            })
            continue

        valid_mask = counts >= min_patients
        valid_means = means[valid_mask]
        rank = int((valid_means > gene_mean).sum() + 1)
        n_valid = int(valid_mask.sum())
        percentile = round(rank / n_valid * 100, 1)

        records.append({
            "dataset": dataset_name,
            "mean_log2FC": gene_mean,
            "rank": rank,
            "n_valid_genes": n_valid,
            "percentile": percentile,
            "n_patients": gene_n_patients,
        })

    return pd.DataFrame(records), real_name


# =========================================================
# 4.3 基因对相关性分析
# =========================================================
def gene_pair_correlation(df: pd.DataFrame, col_meta: pd.DataFrame, gene_a: str, gene_b: str,
                            mode: str = "combined"):
    """计算两个基因在患者层面的log2FC相关性(Pearson r + p值)"""
    match_a = df[df[GENE_COL].astype(str).str.upper() == gene_a.strip().upper()]
    match_b = df[df[GENE_COL].astype(str).str.upper() == gene_b.strip().upper()]

    if len(match_a) == 0:
        return {"found": False, "reason": f"Gene '{gene_a}' not found"}
    if len(match_b) == 0:
        return {"found": False, "reason": f"Gene '{gene_b}' not found"}

    row_a = match_a.iloc[0]
    row_b = match_b.iloc[0]
    real_a = row_a[GENE_COL]
    real_b = row_b[GENE_COL]

    if mode == "combined":
        cols_meta_f = col_meta.copy()
    else:
        cols_meta_f = col_meta[col_meta["severity"] == mode].copy()

    cols = cols_meta_f["column"].tolist()

    paired = pd.DataFrame({
        "column": cols,
        "gene_a_val": pd.to_numeric(row_a[cols].values, errors="coerce"),
        "gene_b_val": pd.to_numeric(row_b[cols].values, errors="coerce"),
    })
    paired = paired.merge(cols_meta_f[["column", "dataset", "severity"]], on="column", how="left")
    paired = paired.dropna(subset=["gene_a_val", "gene_b_val"])

    if len(paired) < 3:
        return {"found": True, "sufficient_data": False, "n_pairs": len(paired),
                "gene_a": real_a, "gene_b": real_b}

    r, p = pearsonr(paired["gene_a_val"], paired["gene_b_val"])

    return {
        "found": True, "sufficient_data": True,
        "gene_a": real_a, "gene_b": real_b,
        "r": r, "p": p, "n_pairs": len(paired),
        "paired_df": paired,
    }


def plot_gene_correlation(corr_result: dict):
    """绘制两个基因log2FC的相关性散点图 + 回归线"""
    if not corr_result["found"] or not corr_result.get("sufficient_data", True):
        return None

    paired = corr_result["paired_df"]
    gene_a = corr_result["gene_a"]
    gene_b = corr_result["gene_b"]
    r = corr_result["r"]
    p = corr_result["p"]
    n = corr_result["n_pairs"]

    fig = px.scatter(
        paired, x="gene_a_val", y="gene_b_val", color="dataset",
        labels={"gene_a_val": f"{gene_a} log2FC", "gene_b_val": f"{gene_b} log2FC"},
        title=f"{gene_a} vs {gene_b}<br><sup>Pearson r = {r:.3f}, p = {p:.2e}, n = {n} patients</sup>",
        template="plotly_white",
    )

    x_vals = paired["gene_a_val"].values
    y_vals = paired["gene_b_val"].values
    slope, intercept = np.polyfit(x_vals, y_vals, 1)
    x_line = np.linspace(x_vals.min(), x_vals.max(), 50)
    y_line = slope * x_line + intercept

    fig.add_trace(go.Scatter(
        x=x_line, y=y_line, mode="lines",
        line=dict(color="black", dash="dash", width=2),
        name="Regression line",
    ))

    fig.add_hline(y=0, line_color="gray", opacity=0.3)
    fig.add_vline(x=0, line_color="gray", opacity=0.3)

    fig.update_layout(height=550)
    return fig


# =========================================================
# 4.4 森林图 (替代 Bar Chart)
# =========================================================
def compute_forest_stats(group_summary: pd.DataFrame):
    """输入group_summary(含group/n/mean/std列)，计算95%CI + Overall汇总行(仿论文做法)"""
    result = group_summary.copy()
    result["se"] = result["std"] / np.sqrt(result["n"])
    result["ci_lower"] = result["mean"] - 1.96 * result["se"]
    result["ci_upper"] = result["mean"] + 1.96 * result["se"]
    result["direction"] = np.where(result["mean"] > 0, "Up", "Down")

    dataset_means = result["mean"].values
    overall_mean = dataset_means.mean()
    overall_se = dataset_means.std(ddof=1) / np.sqrt(len(dataset_means))
    overall_row = pd.DataFrame([{
        "group": "Overall", "n": result["n"].sum(),
        "mean": overall_mean, "median": np.nan, "std": dataset_means.std(ddof=1),
        "se": overall_se,
        "ci_lower": overall_mean - 1.96 * overall_se,
        "ci_upper": overall_mean + 1.96 * overall_se,
        "direction": "Up" if overall_mean > 0 else "Down",
    }])

    return pd.concat([result, overall_row], ignore_index=True)


def plot_forest(query_result: dict):
    """森林图：纵轴数据集(含Overall)，横轴log2FC，mean±95%CI"""
    if not query_result["found"] or query_result["total_valid_n"] == 0:
        return None

    gene_name = query_result["gene_symbol"]
    mode = query_result["mode"]
    forest_stats_df = compute_forest_stats(query_result["group_summary"])

    order = [g for g in forest_stats_df["group"] if g != "Overall"]
    order = sorted(order) + ["Overall"]
    forest_stats_df["group"] = pd.Categorical(forest_stats_df["group"], categories=order[::-1], ordered=True)
    df_plot = forest_stats_df.sort_values("group")

    colors = df_plot["direction"].map({"Up": UP_COLOR, "Down": DOWN_COLOR})
    symbols = ["diamond" if g == "Overall" else "circle" for g in df_plot["group"]]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df_plot["mean"], y=df_plot["group"],
        mode="markers",
        marker=dict(color=colors, size=12, symbol=symbols),
        error_x=dict(
            type="data", symmetric=False,
            array=df_plot["ci_upper"] - df_plot["mean"],
            arrayminus=df_plot["mean"] - df_plot["ci_lower"],
            color="gray", thickness=1.5,
        ),
        text=[f"n={int(n)}" for n in df_plot["n"]],
        hovertemplate="<b>%{y}</b><br>Mean log2FC=%{x:.3f}<br>%{text}<extra></extra>",
    ))

    fig.add_vline(x=0, line_dash="dash", line_color="gray", opacity=0.6)

    mode_display = {"combined": "Severe Dengue (DHF+DSS)", "DHF": "DHF", "DSS": "DSS"}[mode]
    fig.update_layout(
        title=f"{gene_name} — Forest Plot ({mode_display} vs DF)",
        xaxis_title="log2 Fold Change (mean ± 95% CI)",
        yaxis_title="Dataset",
        template="plotly_white",
        height=450,
    )
    return fig
