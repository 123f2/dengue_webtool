import streamlit as st
from utils import load_and_preprocess_data, query_gene, plot_gene_boxplot, plot_forest, GENE_COL

st.set_page_config(page_title="Gene Query — Severe Dengue Explorer", page_icon="🔍", layout="wide")

DATA_PATH = "data/Severe_dengue_DHF_DSS_Xinyang.xlsx"

st.title("🔍 Gene Query")
st.markdown(
    """
    Search for a gene to see how it changes in severe dengue patients compared to mild
    dengue (DF), across 7 independent datasets. All genes are searchable, even those only
    detected in a subset of patients — the number of patients actually analysed is always
    shown alongside the results.
    """
)

try:
    df, col_meta, dataset_summary = load_and_preprocess_data(DATA_PATH)
except FileNotFoundError:
    st.error(f"⚠️ Data file not found at `{DATA_PATH}`. Please add the dataset to the `data/` folder.")
    st.stop()

col_search, col_mode = st.columns([2, 1])

with col_search:
    gene_query = st.text_input(
        "Enter Gene Symbol",
        placeholder="e.g. STAT1 / IFIT3 / CXCL10 / ISG15 / OAS1",
    )

with col_mode:
    analysis_mode = st.selectbox(
        "Select Analysis Group",
        options=["combined", "DHF", "DSS"],
        format_func=lambda m: {
            "combined": "Combined Severe Dengue (DHF + DSS)",
            "DHF": "DHF vs DF",
            "DSS": "DSS vs DF",
        }[m],
    )

st.caption(
    "💡 Example genes to try: **STAT1** (shows cross-dataset heterogeneity), "
    "**IFIT3** (illustrates missing-value handling — not detected in all patients)"
)

st.markdown("---")

if gene_query.strip() == "":
    st.info("👆 Enter a gene symbol above to begin.")
    st.stop()

result = query_gene(df, col_meta, gene_query, mode=analysis_mode)

# if not result["found"]:
#    st.warning(f"Gene **{gene_query}** was not found in this dataset. Please check the spelling.")
#    st.stop()
from difflib import get_close_matches

if not result["found"]:
    st.warning(f"Gene **{gene_query}** was not found in this dataset. Please check the spelling.")
    
    # 新增:模糊匹配提示
    all_gene_names = df[GENE_COL].astype(str).tolist()
    close_matches = get_close_matches(gene_query.upper(), 
                                        [g.upper() for g in all_gene_names], 
                                        n=5, cutoff=0.6)
    if close_matches:
        # 找回原始大小写的基因名(因为上面比对用的是全大写)
        suggestions = [g for g in all_gene_names if g.upper() in close_matches]
        st.info(f"💡 Did you mean: {', '.join(f'**{s}**' for s in suggestions)}?")
    
    st.stop()

if result["total_valid_n"] == 0:
    st.warning(
        f"Gene **{result['gene_symbol']}** was found, but has no detected values under the "
        f"**{analysis_mode}** analysis group. Try a different analysis group."
    )
    st.stop()

# ---- STEP A: Textual results first ----
st.success(f"✅ Found **{result['gene_symbol']}**")
st.subheader("📊 Summary")

detection_rate = result["total_valid_n"] / result["total_n"] * 100
c1, c2, c3 = st.columns(3)
c1.metric("Patients analysed", f"{result['total_valid_n']} / {result['total_n']}")
c2.metric("Detection rate", f"{detection_rate:.0f}%")
c3.metric("Datasets contributing data", result["group_summary"]["group"].nunique())

st.markdown(f"**Analysis group:** {analysis_mode} vs DF")

st.subheader("🧬 Dataset-level Summary Table")
display_table = result["group_summary"].rename(columns={
    "group": "Dataset", "n": "N patients", "mean": "Mean log2FC",
    "median": "Median log2FC", "std": "Std Dev",
})
st.dataframe(
    display_table.style.format({
        "Mean log2FC": "{:.3f}", "Median log2FC": "{:.3f}", "Std Dev": "{:.3f}",
    }),
    use_container_width=True, hide_index=True,
)

# 新增:下载按钮
csv_data = display_table.to_csv(index=False).encode("utf-8")
st.download_button(
    label="⬇️ Download this table as CSV",
    data=csv_data,
    file_name=f"{result['gene_symbol']}_{analysis_mode}_summary.csv",
    mime="text/csv",
)


if analysis_mode != "combined":
    all_datasets = set(col_meta["dataset"].unique())
    covered_datasets = set(result["group_summary"]["group"].unique())
    missing = sorted(all_datasets - covered_datasets)
    if missing:
        st.caption(
            f"ℹ️ Note: {len(missing)} dataset(s) have no {analysis_mode} patients and are "
            f"therefore not shown: {', '.join(missing)}"
        )

st.markdown("---")

# ---- STEP B: Visualization — user chooses chart type ----
st.subheader("📈 Visualisation")

viz_type = st.selectbox(
    "Choose visualisation type",
    options=["Box Plot (all patient data points)", "Forest Plot (mean ± 95% CI)"],
)

if viz_type.startswith("Box Plot"):
    fig = plot_gene_boxplot(result)
    st.caption(
        "Each point is one patient's log2FC value. Box shows median and interquartile range."
    )
else:
    fig = plot_forest(result)
    st.caption(
        "Each point is the mean log2FC for a dataset, with horizontal whiskers showing the "
        "95% confidence interval. The diamond at the bottom (Overall) summarises the mean "
        "log2FC across all datasets, with its own 95% CI based on between-dataset variability — "
        "if this interval does not cross zero, the overall effect is statistically significant."
    )

if fig is not None:
    st.plotly_chart(fig, use_container_width=True)

    # 新增：图表下载按钮
    chart_type_label = "boxplot" if viz_type.startswith("Box Plot") else "forestplot"
    try:
        img_bytes = fig.to_image(format="png", scale=2, width=1000, height=600)
        st.download_button(
            label="⬇️ Download this chart as PNG",
            data=img_bytes,
            file_name=f"{result['gene_symbol']}_{analysis_mode}_{chart_type_label}.png",
            mime="image/png",
        )
    except Exception as e:
        st.caption(
            "ℹ️ Chart image export is temporarily unavailable. "
            "You can still right-click the chart above to save it manually."
        )

else:
    st.warning("No data available to plot for this gene under the current analysis group.")
