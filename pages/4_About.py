import streamlit as st

st.set_page_config(page_title="About — Severe Dengue Explorer", page_icon="ℹ️", layout="wide")

st.title("ℹ️ About This Project")

st.markdown(
    """
    ## Aim

    Severe dengue (DHF and DSS) causes substantial morbidity and mortality, but the host
    gene expression changes behind disease severity are not fully understood. This tool
    helps researchers query and visualize gene expression changes associated with severe
    dengue across several independent patient cohorts.

    ## Data

    Gene expression data (log2 fold-change, severe dengue vs. mild dengue) were collated
    from **7 publicly available datasets**, with **127 severe dengue patients** in total
    (77 DHF, 50 DSS).

    ## Methods

    - **Missing values:** log2FC values are read per patient and per dataset. Missing
      values are not imputed or treated as zero. They are left out of each calculation,
      and every result reports the number of patients analyzed.
    - **Gene Query:** summary statistics are calculated separately for each dataset and
      are **not** pooled across datasets. Results are shown as a summary table, a box
      plot (patient level) and a forest plot (dataset level).
    - **Forest plot:** each dataset's mean log2FC is plotted with its 95% confidence
      interval (mean ± 1.96 × standard error). The **Overall** estimate is the mean of
      the dataset means, with a 95% CI based on the variation between datasets
      (following Chan et al., 2023). If this interval does not cross zero, the gene
      tends to change in one direction across the available datasets.

    ## Gene Rankings & Comparisons

    - **Top Genes per Dataset:** in the chosen dataset, genes detected in at least 3
      patients are ranked by mean log2FC, and the top and bottom N genes are highlighted
      against all other genes. The page can also count how many datasets share the same
      top N genes.
    - **Gene Rank Finder:** shows the mean log2FC, rank and percentile of a gene in each
      dataset (0% = most upregulated, 100% = most downregulated). Percentile is used to
      compare datasets because they rank different numbers of genes. A gene detected in
      fewer than 3 patients in a dataset is reported as "Not detected".
    - **Gene-Gene Correlation:** pairs the log2FC values of two genes in each patient
      (pooled across the datasets in the selected group) and gives Pearson's r and the
      p-value. Because log2FC is relative to mild dengue, this shows whether two genes
      change together among severe dengue patients, not whether they are associated with
      severity.

    ## Limitations

    - Sample sizes are small: 10 to 32 patients per dataset, and as few as 6 in some DHF
      or DSS subgroups.
    - Datasets differ in platform, patient population and number of genes (GSE25226
      covers only about 1,650 genes).
    - Sex-linked genes (e.g. *XIST*, *RPS4Y1*, *UTY*) appear among the top-ranked genes
      in some datasets. Sex information is not included in this tool, so such genes may
      reflect differences in sex composition rather than disease severity.
    - Ranks and percentiles are exploratory: they are computed separately for each
      dataset and are not adjusted for multiple testing.
    - Gene-gene correlations pool patients across datasets without batch correction, so
      a strong correlation can be partly driven by differences between datasets.
    """
)

st.markdown("---")
st.page_link("home.py", label="← Back to Home", icon="🏠")
