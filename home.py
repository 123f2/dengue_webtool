import streamlit as st
from pathlib import Path
from utils import load_and_preprocess_data

st.set_page_config(
    page_title="Severe Dengue Transcriptomic Explorer",
    page_icon="🧬",
    layout="wide",
)

DATA_PATH = "data/Severe_dengue_DHF_DSS_Xinyang.xlsx"

# ============ 新增：顶部banner图片 ============
BANNER_PATH = "assets/dengue_banner.jpg"
banner_path_obj = Path(BANNER_PATH)

if banner_path_obj.exists():
    col_left, col_center, col_right = st.columns([1, 6, 1])
    with col_center:
        st.image(str(banner_path_obj), use_container_width=True)
    st.caption(
        "<div style='text-align: center; color: #8a99a8; font-size: 0.85rem; margin-top: 0.3rem;'>"
        "Sources: WHO Dengue fact sheet (21 Aug 2025); WHO Global dengue surveillance dashboard; "
        "NEA Dengue Cases page (accessed 6 Oct 2026)</div>",
        unsafe_allow_html=True,
    )
    st.write("")
else:
    st.info("ℹ️ Banner image not found — place it at `assets/dengue_banner.jpg` to display it here.")
# ============ 插入结束 ============


st.markdown(
    """
    <style>
    .main-header { font-size: 2.3rem; font-weight: 700; color: #1a3c5e; margin-bottom: 0rem; }
    .sub-header { font-size: 1.1rem; color: #5a6c7d; margin-top: 0rem; margin-bottom: 1.5rem; }
    
    /* 新增：给banner图片加圆角+阴影 */
    div[data-testid="stImage"] img {
        border-radius: 12px;
        box-shadow: 0 4px 16px rgba(0, 0, 0, 0.12);
    }
    
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

st.markdown(
    """
    Dengue is a mosquito-borne viral infection that is now endemic in more than 100 countries.
    WHO recorded 14.6 million cases in 2024, the highest total on record, and its global
    surveillance dashboard lists 7.86 million cases between January 2025 and August 2026.
    In Singapore, the National Environment Agency reported 2,564 cases in the first 39
    epidemiological weeks of 2026. Most infections are mild, but a minority develop severe
    dengue, for which there is no specific treatment.
    """
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

col3, col4 = st.columns(2)
with col3:
    st.subheader("📊 Gene Rankings")
    st.write("Find top genes per dataset, look up a gene's rank across datasets, or correlate two genes.")
    st.page_link("pages/3_Gene_Rankings.py", label="Open Gene Rankings", icon="📊")

with col4:
    st.subheader("ℹ️ About")
    st.write("Learn about the project's aim, data sources, methodology, and limitations.")
    st.page_link("pages/4_About.py", label="Open About", icon="ℹ️")

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
