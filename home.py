import streamlit as st
from utils import load_and_preprocess_data

st.set_page_config(
    page_title="Severe Dengue Transcriptomic Explorer",
    page_icon="🧬",
    layout="wide",
)

DATA_PATH = "data/Severe_dengue_DHF_DSS_Xinyang.xlsx"

st.markdown(
    """
    <style>
    .main-header { font-size: 2.3rem; font-weight: 700; color: #1a3c5e; margin-bottom: 0rem; }
    .sub-header { font-size: 1.1rem; color: #5a6c7d; margin-top: 0rem; margin-bottom: 1.5rem; }
    </style>
    """,
    unsafe_allow_html=True,
)

st.markdown('<p class="main-header">🧬 Severe Dengue Transcriptomic Explorer</p>', unsafe_allow_html=True)
st.markdown(
    '<p class="sub-header">A resource for querying gene expression changes in severe dengue '
    '(DHF / DSS) relative to mild dengue (DF), integrated across 7 published patient cohorts.</p>',
    unsafe_allow_html=True,
)

st.markdown("---")

try:
    df, col_meta, dataset_summary = load_and_preprocess_data(DATA_PATH)
    data_loaded = True
except FileNotFoundError:
    data_loaded = False
    st.error(f"⚠️ Data file not found at `{DATA_PATH}`. Please place the dataset in the `data/` folder.")

if data_loaded:
    n_genes = df.shape[0]
    n_patients = col_meta.shape[0]
    n_datasets = col_meta["dataset"].nunique()
    n_dhf = int((col_meta["severity"] == "DHF").sum())
    n_dss = int((col_meta["severity"] == "DSS").sum())

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Genes covered", f"{n_genes:,}")
    c2.metric("Patients", n_patients)
    c3.metric("Independent datasets", n_datasets)
    c4.metric("DHF / DSS split", f"{n_dhf} / {n_dss}")

st.markdown("---")

st.markdown(
    """
    ## What this tool does

    This webtool lets you search for **any gene** and see how its expression changes in
    severe dengue patients (DHF and/or DSS) compared to mild dengue (DF), across 7
    independently published patient cohorts.

    - 🔍 **Query any gene** — results shown even if the gene was only detected in a subset
      of patients (no genes are excluded).
    - 📊 **See the real distribution** — patient-level log2 fold-change values are shown
      per dataset, not just a single averaged number.
    - 🧮 **Choose your comparison** — combined severe dengue (DHF+DSS), DHF-only, or
      DSS-only.

    ## How the data is organised

    Gene expression is provided as **log2 fold-change (log2FC)** values, already
    calculated per patient within each original study. Values are integrated — but
    **not merged or averaged** — across 7 datasets, so dataset-level patterns remain
    visible.

    Empty cells mean the gene was **not detected** in that patient/platform — treated
    as missing data (not zero), excluded only from that specific calculation. The
    number of patients actually analysed is always shown alongside results.
    """
)

if data_loaded:
    with st.expander("📋 View dataset composition (patients per study)"):
        st.dataframe(dataset_summary, use_container_width=True, hide_index=True)

# ============ 在这里插入 Quick Search 代码 ============
st.markdown("---")
st.markdown("## Quick Search")

col_q1, col_q2 = st.columns([3, 1])
with col_q1:
    quick_gene = st.text_input(
        "Search a gene now",
        placeholder="e.g. STAT1",
        label_visibility="collapsed",
        key="home_quick_search",
    )
with col_q2:
    quick_search_clicked = st.button("🔍 Search", use_container_width=True)

if quick_search_clicked:
    st.switch_page("pages/1_Gene_Query.py")
# ============ 插入结束 ============


st.markdown("---")
st.markdown("## Get started")

col1, col2 = st.columns(2)
with col1:
    st.subheader("🔍 Gene Query")
    st.write("Search for a gene and view its log2FC distribution across all 7 datasets.")
    st.page_link("pages/1_Gene_Query.py", label="Open Gene Query", icon="🔍")

with col2:
    st.subheader("🧪 Dataset Explorer")
    st.write("Browse the 7 underlying studies and their patient composition.")
    st.page_link("pages/2_Dataset_Explorer.py", label="Open Dataset Explorer", icon="🧪")

st.markdown("---")

with st.expander("🛣️ Roadmap — coming in future updates"):
    st.markdown(
        """
        - Automated ranking of the most consistently differentially expressed genes
          across datasets (candidate drug targets)
        - Stratification by patient demographics (age, sex) where available
        - Pathway / gene-set level summaries
        """
    )

st.markdown("---")
st.caption("Built with Streamlit · Data integrated from 7 publicly available dengue transcriptomic datasets")
