import streamlit as st
from utils import load_and_preprocess_data

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