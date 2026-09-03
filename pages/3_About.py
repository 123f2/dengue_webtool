import streamlit as st

st.set_page_config(page_title="About — Severe Dengue Explorer", page_icon="ℹ️", layout="wide")

st.title("ℹ️ About This Project")

st.markdown(
    """
    ## Aim

    Severe dengue (DHF and DSS) causes substantial morbidity and mortality, and the host
    transcriptomic changes underlying disease severity are not fully understood. This tool
    was built to help researchers quickly query and visualise gene expression changes
    associated with severe dengue, integrated across multiple independently published
    patient cohorts.

    ## Data

    Gene expression data (log2 fold-change, severe dengue vs. mild dengue) were collated
    from **7 publicly available datasets**, comprising **127 severe dengue patients**
    (77 DHF, 50 DSS) in total.

    ## How missing values are handled

    Not every gene was detected on every microarray platform or in every patient. Rather
    than imputing missing values or treating them as zero, this tool **excludes only the
    specific missing data points** from each calculation, and always reports how many
    patients were actually analysed.

    ## Methodology summary

    1. Raw log2FC values are read per patient, per dataset.
    2. For a queried gene, all non-missing values are retrieved and grouped by dataset
       (and, if selected, filtered to DHF-only or DSS-only comparisons).
    3. Dataset-level summary statistics are calculated independently for each dataset —
       values are **not** pooled across datasets into a single combined sample.
    4. Results are shown as both a summary table and an interactive plot chosen by
       the user.

    ## Limitations

    - Sample sizes per dataset are modest (as low as 6 patients in some subgroups).
    - Datasets differ in platform, patient population, and demographic composition.
    - This version does not yet correct for demographic differences or provide formal
      statistical significance testing.

    ## Contact

    For questions about the underlying data or methodology, please reach out to the
    project team.
    """
)