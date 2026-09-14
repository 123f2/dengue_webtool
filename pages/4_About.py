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

    ## Gene Rankings & Comparisons (advanced features)

    Beyond single-gene queries, the **Gene Rankings** page offers three additional ways
    to explore the data:

    - **Top Genes per Dataset** — for a chosen dataset, all genes with at least 3
      detected patients are ranked by mean log2FC. The top and bottom N genes (N is
      user-selectable) are highlighted against the full distribution of genes, so
      extreme values can be seen in context rather than in isolation.
    - **Gene Rank Finder** — for a gene of interest, this shows its rank, percentile,
      and mean log2FC within each dataset. Percentile (0% = most upregulated, 100% =
      most downregulated) is the metric intended for comparing a gene's relative
      standing *across* datasets, since datasets differ in how many genes pass the
      minimum-detection threshold and therefore have different total gene counts.
    - **Gene-Gene Correlation** — for any two genes, their log2FC values are paired at
      the patient level (using patients from the selected analysis group who have data
      for both genes) and correlated using Pearson's r, with the associated p-value.
      This tests whether two genes tend to move together across patients, which can
      support (but does not on its own prove) a shared biological pathway or mechanism.

    ## Forest plots

    In place of a bar chart, gene-level results can also be shown as a **forest plot**:
    each dataset's mean log2FC is plotted with its 95% confidence interval (mean ±
    1.96 × standard error). An additional **Overall** estimate is shown at the bottom,
    calculated as the mean of the per-dataset means, with its own 95% CI reflecting the
    variability *between* datasets — if this interval does not cross zero, the gene's
    overall direction of change is considered statistically supported across the
    available datasets. This approach follows standard meta-analysis forest plot
    conventions and avoids the visual distortion of bar charts, which must start at zero.

    ## Limitations

    - Sample sizes per dataset are modest (as low as 6 patients in some subgroups).
    - Datasets differ in platform, patient population, and demographic composition —
      in particular, some datasets show an imbalance in sex ratio between severity
      groups, and sex-linked genes (e.g. *XIST*, *RPS4Y1*, *UTY*) have been observed
      among the top-ranked genes in some datasets. This is a reminder that a gene's
      prominence in this tool does not by itself establish that it is driven by disease
      severity rather than cohort composition.
    - Ranks and percentiles are computed independently per dataset and are not adjusted
      for multiple testing; they are intended as an exploratory guide rather than a
      confirmatory statistical result.
    - Gene-gene correlations are computed by pooling patients across datasets within
      the selected analysis group; they do not account for between-dataset batch
      effects, so a strong correlation could in principle be partly driven by
      between-dataset differences rather than a purely within-patient relationship.
    - This version does not yet correct for demographic differences or provide formal
      multiple-testing correction.

    ## Contact

    For questions about the underlying data or methodology, please reach out to the
    project team.
    """
)