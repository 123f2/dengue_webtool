import streamlit as st
from utils import (
    load_and_preprocess_data,
    get_overall_dataset_distribution,
    plot_overall_distribution,
    GENE_COL,
)

st.set_page_config(page_title="Dataset Explorer — Severe Dengue Explorer", page_icon="🧪", layout="wide")

DATA_PATH = "data/Severe_dengue_DHF_DSS_Xinyang.xlsx"

st.title("🧪 Dataset Explorer")
st.markdown(
    """
    This tool integrates **7 independently published dengue transcriptomic datasets**.
    Each dataset contributes patient-level log2 fold-change (log2FC) values comparing
    severe dengue (DHF and/or DSS) against mild dengue (DF). Datasets are kept separate
    throughout the analysis — values are never pooled into one combined sample.
    """
)

try:
    df, col_meta, dataset_summary = load_and_preprocess_data(DATA_PATH)
except FileNotFoundError:
    st.error(f"⚠️ Data file not found at `{DATA_PATH}`. Please add the dataset to the `data/` folder.")
    st.stop()

# =========================
# 1. 数据集患者组成
# =========================
st.subheader("📋 Patients per dataset")

pivot = dataset_summary.pivot(index="dataset", columns="severity", values="n_patients").fillna(0).astype(int)
pivot["Total"] = pivot.sum(axis=1)
pivot = pivot.reset_index().rename(columns={"dataset": "Dataset"})

st.dataframe(pivot, use_container_width=True, hide_index=True)

c1, c2, c3 = st.columns(3)
c1.metric("Total datasets", col_meta["dataset"].nunique())
c2.metric("Total DHF patients", int((col_meta["severity"] == "DHF").sum()))
c3.metric("Total DSS patients", int((col_meta["severity"] == "DSS").sum()))

st.markdown("---")

# =========================
# 2. 原始数据预览
# =========================
st.subheader("📄 Raw Data Preview")
st.caption(
    "A preview of the underlying expression matrix. Rows are genes; columns are "
    "individual patient comparisons (log2FC vs. DF). Empty cells indicate the gene "
    "was not detected in that patient/platform."
)

all_datasets = sorted(col_meta["dataset"].unique())
preview_dataset = st.selectbox("Select a dataset to preview its columns", options=all_datasets)

preview_cols = col_meta[col_meta["dataset"] == preview_dataset]["column"].tolist()
n_show_rows = st.slider("Number of gene rows to preview", min_value=5, max_value=50, value=10, step=5)

preview_df = df[[GENE_COL] + preview_cols].head(n_show_rows)
st.dataframe(
    preview_df,
    use_container_width=True,
    hide_index=True,
    height=300,
)

st.caption(f"Showing {n_show_rows} of {df.shape[0]:,} total genes, "
           f"{len(preview_cols)} patient column(s) from {preview_dataset}.")


# =========================
# 4. 备注
# =========================
st.subheader("ℹ️ Notes on dataset composition")
st.markdown(
    """
    - **GSE43777** in this table combines two platform-specific column groups (GPL201 and
      GPL570) from the same original study, which profiled patients using two different
      microarray generations.
    - Not all datasets include both DHF and DSS comparisons — some studies only reported
      one severity subtype relative to DF. When filtering to "DHF only" or "DSS only" on
      the Gene Query page, datasets without that comparison type are automatically excluded.
    - Patient demographic details (age, sex distribution) vary across datasets and are
      documented in the original source publications. Some cohorts show demographic
      imbalance between severity groups, which is worth keeping in mind when interpreting
      cross-cohort consistency.
    """
)

st.markdown("---")
st.page_link("home.py", label="← Back to Home", icon="🏠")
