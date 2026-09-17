import streamlit as st
import pandas as pd
from utils import (
    load_and_preprocess_data,
    get_gene_landscape,
    mark_top_hits,
    plot_gene_rank_scatter,
    get_top_hits_table,
    gene_rank_across_datasets,
    gene_pair_correlation,
    plot_gene_correlation,
)

st.set_page_config(page_title="Gene Rankings — Severe Dengue Explorer", page_icon="📊", layout="wide")

DATA_PATH = "data/Severe_dengue_DHF_DSS_Xinyang.xlsx"

st.title("📊 Gene Rankings & Comparisons")
st.markdown(
    """
    Explore which genes stand out the most within a specific dataset, look up how a gene
    of interest ranks across all datasets, or compare two genes to see if their expression
    changes are related.
    """
)

try:
    df, col_meta, dataset_summary = load_and_preprocess_data(DATA_PATH)
except FileNotFoundError:
    st.error(f"⚠️ Data file not found at `{DATA_PATH}`. Please add the dataset to the `data/` folder.")
    st.stop()

MODE_LABELS = {
    "combined": "Combined Severe Dengue (DHF + DSS)",
    "DHF": "DHF vs DF",
    "DSS": "DSS vs DF",
}

tab1, tab2, tab3 = st.tabs(["🎯 Top Genes per Dataset", "🔍 Gene Rank Finder", "🔗 Gene-Gene Correlation"])

# =========================================================
# TAB 1: 数据集内Top基因排名 (4.1)
# =========================================================
with tab1:
    st.subheader("Top upregulated / downregulated genes in a specific dataset")

    c1, c2, c3 = st.columns([2, 1, 1])
    with c1:
        all_datasets = sorted(col_meta["dataset"].unique())
        selected_dataset = st.selectbox("Select a dataset", options=all_datasets, key="t1_dataset")
    with c2:
        selected_mode_t1 = st.selectbox(
            "Analysis Group", options=["combined", "DHF", "DSS"],
            format_func=lambda m: MODE_LABELS[m], key="t1_mode",
        )
    with c3:
        top_n = st.slider("Number of top hits", min_value=1, max_value=50, value=10, key="t1_topn")

    landscape = get_gene_landscape(df, col_meta, selected_dataset, mode=selected_mode_t1)

    if landscape is None:
        st.warning(
            f"⚠️ {selected_dataset} has no {selected_mode_t1} patients — "
            f"please choose a different analysis group."
        )
    else:
        marked = mark_top_hits(landscape, top_n=top_n)
        fig = plot_gene_rank_scatter(marked, selected_dataset, mode_label=MODE_LABELS[selected_mode_t1])
        st.plotly_chart(fig, use_container_width=True)
        st.caption(
            "Gray points are all other genes; hover over a highlighted point to see the gene name. "
            f"{len(landscape):,} genes were included in this ranking "
            "(genes detected in fewer than 3 patients within this dataset are excluded)."
        )

        up_table, down_table = get_top_hits_table(marked)

        col_up, col_down = st.columns(2)
        with col_up:
            st.markdown(f"**Top {top_n} Upregulated**")
            st.dataframe(
                up_table.style.format({"Mean log2FC": "{:.3f}"}),
                use_container_width=True, hide_index=True, height=min(35 * (top_n + 1), 500),
            )
            csv_up = up_table.to_csv(index=False).encode("utf-8")
            st.download_button(
                "⬇️ Download Top Upregulated as CSV", data=csv_up,
                file_name=f"{selected_dataset}_{selected_mode_t1}_top{top_n}_up.csv", mime="text/csv",
                key="t1_dl_up",
            )
        with col_down:
            st.markdown(f"**Top {top_n} Downregulated**")
            st.dataframe(
                down_table.style.format({"Mean log2FC": "{:.3f}"}),
                use_container_width=True, hide_index=True, height=min(35 * (top_n + 1), 500),
            )
            csv_down = down_table.to_csv(index=False).encode("utf-8")
            st.download_button(
                "⬇️ Download Top Downregulated as CSV", data=csv_down,
                file_name=f"{selected_dataset}_{selected_mode_t1}_top{top_n}_down.csv", mime="text/csv",
                key="t1_dl_down",
            )

# =========================================================
# TAB 2: 基因跨数据集排名查询 (4.2)
# =========================================================
with tab2:
    st.subheader("Look up a gene's rank and log2FC across all datasets")

    c1, c2 = st.columns([2, 1])
    with c1:
        gene_query_t2 = st.text_input(
            "Enter Gene Symbol", placeholder="e.g. STAT1 / UQCRQ", key="t2_gene",
        )
    with c2:
        selected_mode_t2 = st.selectbox(
            "Analysis Group", options=["combined", "DHF", "DSS"],
            format_func=lambda m: MODE_LABELS[m], key="t2_mode",
        )

    if gene_query_t2.strip():
        rank_df, real_name = gene_rank_across_datasets(df, col_meta, gene_query_t2, mode=selected_mode_t2)

        if rank_df is None:
            st.warning(f"Gene **{gene_query_t2}** was not found. Please check the spelling.")
        else:
            st.success(f"✅ Found **{real_name}**")

            display_df = rank_df.rename(columns={
                "dataset": "Dataset", "mean_log2FC": "Mean log2FC", "rank": "Rank",
                "n_valid_genes": "Total Genes Ranked", "percentile": "Percentile",
                "n_patients": "N Patients",
            }).copy()
            display_df["Percentile"] = display_df["Percentile"].apply(
                lambda p: f"{p}%" if pd.notna(p) else "—"
            )
            display_df["Mean log2FC"] = display_df["Mean log2FC"].apply(
                lambda v: f"{v:.3f}" if pd.notna(v) else "Not detected"
            )
            display_df["Rank"] = display_df["Rank"].apply(
                lambda r: f"{int(r):,}" if pd.notna(r) else "—"
            )
            display_df["Total Genes Ranked"] = display_df["Total Genes Ranked"].apply(
               lambda v: f"{int(v):,}" if pd.notna(v) else "—"
            )

            st.dataframe(display_df, use_container_width=True, hide_index=True)
            st.caption(
                "Percentile: 0% = most upregulated gene in that dataset, 100% = most downregulated. "
                "Datasets where the gene was detected in fewer than 3 patients show 'Not detected'."
            )

            csv_data = display_df.to_csv(index=False).encode("utf-8")
            st.download_button(
                "⬇️ Download this table as CSV", data=csv_data,
                file_name=f"{real_name}_{selected_mode_t2}_rank_across_datasets.csv", mime="text/csv",
                key="t2_dl",
            )
    else:
        st.info("👆 Enter a gene symbol above to see its rank across all datasets.")

# =========================================================
# TAB 3: 基因对相关性分析 (4.3)
# =========================================================
with tab3:
    st.subheader("Correlate two genes' log2FC values across patients")
    st.caption(
        "Select two genes to see whether their expression changes tend to move together "
        "across patients (pooled from the selected analysis group)."
    )

    c1, c2, c3 = st.columns([2, 2, 1])
    with c1:
        gene_a_input = st.text_input("Gene A", placeholder="e.g. STAT1", key="t3_gene_a")
    with c2:
        gene_b_input = st.text_input("Gene B", placeholder="e.g. IFIT3", key="t3_gene_b")
    with c3:
        selected_mode_t3 = st.selectbox(
            "Analysis Group", options=["combined", "DHF", "DSS"],
            format_func=lambda m: MODE_LABELS[m], key="t3_mode",
        )

    if gene_a_input.strip() and gene_b_input.strip():
        corr_result = gene_pair_correlation(df, col_meta, gene_a_input, gene_b_input, mode=selected_mode_t3)

        if not corr_result["found"]:
            st.warning(f"⚠️ {corr_result['reason']}")
        elif not corr_result.get("sufficient_data", True):
            st.warning(
                f"Found both genes, but only {corr_result['n_pairs']} patients have data for "
                f"both genes under the **{selected_mode_t3}** group — not enough to compute a reliable correlation."
            )
        else:
            fig = plot_gene_correlation(corr_result)
            st.plotly_chart(fig, use_container_width=True)

            r = corr_result["r"]
            p = corr_result["p"]
            n = corr_result["n_pairs"]
            significance = "statistically significant" if p < 0.05 else "not statistically significant"
            st.caption(
                f"Pearson r = {r:.3f}, p = {p:.2e}, based on {n} patients with data for both genes "
                f"— this correlation is **{significance}** (α = 0.05). "
                f"Each point is one patient; colors indicate which dataset the patient came from."
            )
    else:
        st.info("👆 Enter both gene symbols above to see their correlation.")
