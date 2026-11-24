# Code release header: turns the scib-metrics results of one seed into raw and min-max scaled score tables and diagnostic plots.
# Input: scib_results_prefix.csv of one seed (from 03_run_python_integration.py or R1_rescore_scib.py). BBKNN's PCR comparison is set to NaN.
# Batch correction and biological conservation are means of five metrics each (NaN skipped). Total = 0.4 x batch correction + 0.6 x biological conservation.
# Each metric is min-max scaled across methods. A step not included averages seeds 0 and 1. Output names use an internal numbering (Table1-Table4 CSVs and Fig8 plots). Methods: Benchmarking of integration methods.
"""04_scib_tables_figures.py - BBKNN fix, min-max scaling, Tables 1-4 CSVs and Figure 8 panels.

Steps:
  * BBKNN† 'PCR comparison' -> NaN (graph-only method)
  * Bio conservation = mean of 5 bio metrics; Batch correction = mean of 5 batch metrics (NaN skipped);
    Total = 0.4 x Batch + 0.6 x Bio
  * min-max: each of the 10 metric columns scaled to [0,1] across methods (NaN ignored), aggregates recomputed
    (this is NOT scib-metrics' own get_results(min_max_scale=True))
  * Figure 8a/b = scatter of min-max Batch vs Bio; Figure 8c/d = plottable table, numbers min-max, circle color raw

Usage:
  python 04_scib_tables_figures.py --mode allcells|luminal --prefix-csv <tables/scib_results_prefix.csv> --outdir <dir>
         [--compare-to <raw or min-max table csv>]
  python 04_scib_tables_figures.py --mode allcells|luminal --check-published --outdir <dir>
         (re-applies the fix to the source table PUB['src'] and must reproduce the PUB raw and min-max tables)
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
import matplotlib.colors as mcolors
from plottable import Table, ColumnDefinition
from plottable.plots import bar
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

p = argparse.ArgumentParser()
p.add_argument('--mode', choices=['allcells', 'luminal'], required=True)
p.add_argument('--prefix-csv', default=None)
p.add_argument('--outdir', required=True)
p.add_argument('--check-published', action='store_true')
p.add_argument('--compare-to', default=None, help='published raw or min-max table to diff against (same kind)')
p.add_argument('--tag', default='')
a = p.parse_args()

Z = '<BASE_DIR>/17datasets/Zenodo_file_submission/supplementary_tables'
PY = '<BASE_DIR>/17datasets/python_final'
PUB = {'allcells': dict(raw=f'{Z}/Suppl_Table_4_all_raw_scib.csv', mm=f'{Z}/Suppl_Table_5_all_norm_scib.csv',
                        src=f'{PY}/scib_a17_hit1_50pc2000g/scib_results.csv', t_raw='Table1_allcells_raw_scib.csv',
                        t_mm='Table2_allcells_minmax_scib.csv', fig='Fig8a_scatter_minmax_allcells',
                        tab='Fig8c_results_table_minmax_rawcolour_allcells', name='all cell types (50pc2000g)'),
       'luminal': dict(raw=f'{Z}/Suppl_Table_6_lum_raw_scib.csv', mm=f'{Z}/Suppl_Table_7_lum_norm_scib.csv',
                       src=f'{PY}/032025_scib/scib_lum_v1_30pc3000g/scib_results.csv', t_raw='Table3_luminal_raw_scib.csv',
                       t_mm='Table4_luminal_minmax_scib.csv', fig='Fig8b_scatter_minmax_luminal',
                       tab='Fig8d_results_table_minmax_rawcolour_luminal', name='luminal cells (30pc3000g)')}[a.mode]
OUT = a.outdir; os.makedirs(OUT, exist_ok=True)
RAW_CSV = PUB['src'] if a.check_published else a.prefix_csv
TAG = a.tag or ('published tables re-processed (check)' if a.check_published else 're-run')

BBKNN_LABEL = 'BBKNN†'
BBKNN_NA_METRICS = ['PCR comparison']
df_raw = pd.read_csv(RAW_CSV, index_col=0)
metric_type = df_raw.loc['Metric Type']
df_scores = df_raw.drop(index='Metric Type').apply(pd.to_numeric, errors='coerce')
bio_cols = metric_type[metric_type == 'Bio conservation'].index.tolist()
batch_cols = metric_type[metric_type == 'Batch correction'].index.tolist()
agg_cols = metric_type[metric_type == 'Aggregate score'].index.tolist()
indiv_cols = bio_cols + batch_cols
assert len(bio_cols) == 5 and len(batch_cols) == 5, (bio_cols, batch_cols)

df_fixed = df_scores.copy()
for col in BBKNN_NA_METRICS:
    if BBKNN_LABEL in df_fixed.index:
        print(f'{BBKNN_LABEL} [{col}]: {df_fixed.loc[BBKNN_LABEL, col]} -> NaN')
        df_fixed.loc[BBKNN_LABEL, col] = np.nan

def recompute_aggregates(df, bio_cols, batch_cols):
    df = df.copy()
    df['Bio conservation'] = df[bio_cols].mean(axis=1, skipna=True)
    df['Batch correction'] = df[batch_cols].mean(axis=1, skipna=True)
    df['Total'] = 0.4 * df['Batch correction'] + 0.6 * df['Bio conservation']
    return df
df_fixed = recompute_aggregates(df_fixed, bio_cols, batch_cols)

def minmax_scale(df, cols):
    df_s = df.copy()
    for col in cols:
        mn, mx = df[col].min(skipna=True), df[col].max(skipna=True)
        df_s[col] = (df[col] - mn) / (mx - mn) if mx > mn else 0.0
    return df_s
df_scaled = recompute_aggregates(minmax_scale(df_fixed, indiv_cols), bio_cols, batch_cols)

def with_metric_type(df):
    return pd.concat([df, pd.DataFrame([metric_type], index=['Metric Type'])])
with_metric_type(df_fixed).to_csv(os.path.join(OUT, PUB['t_raw']))
with_metric_type(df_scaled).to_csv(os.path.join(OUT, PUB['t_mm']))
print('raw aggregates\n', df_fixed[agg_cols].round(4).sort_values('Total', ascending=False).to_string())
print('min-max aggregates\n', df_scaled[agg_cols].round(4).sort_values('Total', ascending=False).to_string())

if a.check_published:
    for nm, mine, ref in [('raw', df_fixed, PUB['raw']), ('min-max', df_scaled, PUB['mm'])]:
        r = pd.read_csv(ref, index_col=0).drop(index='Metric Type').apply(pd.to_numeric, errors='coerce')
        d = (mine.loc[r.index, r.columns] - r).abs().max().max()
        print(f'CHECK vs {os.path.basename(ref)} ({nm}): max |diff| = {d:.3g}')
if a.compare_to:
    r = pd.read_csv(a.compare_to, index_col=0).drop(index='Metric Type').apply(pd.to_numeric, errors='coerce')
    mine = df_scaled if ('norm' in a.compare_to or 'minmax' in a.compare_to) else df_fixed
    common = [i for i in r.index if i in mine.index]
    delta = (mine.loc[common, r.columns] - r.loc[common]).round(4)
    delta.to_csv(os.path.join(OUT, 'delta_vs_published.csv'))
    print('new - published\n', delta[agg_cols].to_string())

FOOTNOTE = ('† BBKNN: graph-only method — N/A metrics excluded from aggregate (Luecken et al. 2022, Nat Methods).  '
            'Total = 0.4 × Batch + 0.6 × Bio.')
def normed_cmap(series, cmap=mcm.PRGn, num_stds=2.5):
    s = pd.to_numeric(series, errors='coerce').dropna()
    norm = mcolors.Normalize(vmin=s.median() - num_stds * s.std(), vmax=s.median() + num_stds * s.std())
    return mcm.ScalarMappable(norm=norm, cmap=cmap).to_rgba
def save_fig(fig, stem):
    for ext, dpi in (('png', 150), ('pdf', None)):
        kw = dict(bbox_inches='tight')
        if dpi: kw['dpi'] = dpi
        fig.savefig(os.path.join(OUT, f'{stem}.{ext}'), **kw)
    plt.close('all')

def make_results_table(df_data, title, out_stem, color_df=None):
    mt = metric_type.to_dict()
    plot_df = df_data.sort_values(by='Total', ascending=False).astype(object)
    plot_df.insert(0, 'Method', plot_df.index)
    bio_c = [c for c in df_data.columns if mt.get(c) == 'Bio conservation']
    batch_c = [c for c in df_data.columns if mt.get(c) == 'Batch correction']
    agg_c = [c for c in df_data.columns if mt.get(c) == 'Aggregate score']
    df_num = (df_data if color_df is None else color_df)[bio_c + batch_c].apply(pd.to_numeric, errors='coerce')
    col_defs = [ColumnDefinition('Method', width=1.5, textprops={'ha': 'left', 'weight': 'bold'})]
    for col in bio_c + batch_c:
        cmap_fn = normed_cmap(df_num[col])
        group = 'Bio conservation' if col in bio_c else 'Batch correction'
        def _fmt(v):
            try:
                f = float(v); return 'N/A' if np.isnan(f) else f'{f:.2f}'
            except Exception: return 'N/A'
        if color_df is None:
            def _cmap(v, fn=cmap_fn):
                try:
                    f = float(v); return (0.85, 0.85, 0.85, 1.0) if np.isnan(f) else fn(f)
                except Exception: return (0.85, 0.85, 0.85, 1.0)
        else:   # circle color = raw value of the same method (looked up by the displayed scaled value)
            lut = dict(zip([round(float(x), 12) if pd.notna(x) else 'nan' for x in pd.to_numeric(df_data[col], errors='coerce')],
                           pd.to_numeric(color_df[col], errors='coerce')))
            def _cmap(v, fn=cmap_fn, lut=lut):
                try:
                    f = float(v)
                    if np.isnan(f): return (0.85, 0.85, 0.85, 1.0)
                    r = lut.get(round(f, 12), np.nan)
                    return (0.85, 0.85, 0.85, 1.0) if np.isnan(r) else fn(r)
                except Exception: return (0.85, 0.85, 0.85, 1.0)
        col_defs.append(ColumnDefinition(col, title=col.replace(' ', '\n', 1), width=1,
                                         textprops={'ha': 'center', 'bbox': {'boxstyle': 'circle', 'pad': 0.25}, 'fontsize': 8},
                                         cmap=_cmap, group=group, formatter=_fmt))
    for i, col in enumerate(agg_c):
        col_defs.append(ColumnDefinition(col, width=1, title=col.replace(' ', '\n', 1), plot_fn=bar,
                                         plot_kw={'cmap': mcm.YlGnBu, 'plot_bg_bar': False, 'annotate': True, 'height': 0.9,
                                                  'formatter': '{:.2f}'},
                                         group=mt.get(col, 'Aggregate score'), border='left' if i == 0 else None))
    with matplotlib.rc_context({'svg.fonttype': 'none'}):
        fig, ax = plt.subplots(figsize=(len(col_defs) * 1.25, 3 + 0.3 * len(plot_df)))
        Table(plot_df, cell_kw={'linewidth': 0, 'edgecolor': 'k'}, column_definitions=col_defs, ax=ax, row_dividers=True,
              footer_divider=True, textprops={'fontsize': 9}, row_divider_kw={'linewidth': 1, 'linestyle': (0, (1, 5))},
              col_label_divider_kw={'linewidth': 1, 'linestyle': '-'}, column_border_kw={'linewidth': 1, 'linestyle': '-'},
              index_col='Method')
    ax.set_title(title, fontsize=10, pad=12)
    fig.text(0.01, -0.05, FOOTNOTE, fontsize=7.5, color='#555', ha='left')
    save_fig(fig, out_stem)

def make_scatter(df, title, out_stem):
    bc, bio = 'Batch correction', 'Bio conservation'
    df_p = df[[bc, bio]].apply(pd.to_numeric, errors='coerce').dropna()
    colors = plt.cm.tab10.colors
    cmap_m = {m: colors[i % len(colors)] for i, m in enumerate(df.index)}
    fig, ax = plt.subplots(figsize=(6, 5))
    for method, row in df_p.iterrows():
        is_bb = 'BBKNN' in str(method)
        ax.scatter(row[bc], row[bio], s=160 if is_bb else 100, marker='*' if is_bb else 'o', color=cmap_m[method],
                   zorder=3, label=method, edgecolors='k' if is_bb else 'none', linewidths=0.5)
        ax.annotate(method, (row[bc], row[bio]), textcoords='offset points', xytext=(6, 3), fontsize=8)
    ax.set_xlabel('Batch correction', fontsize=12); ax.set_ylabel('Bio conservation', fontsize=12)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axhline(0.5, color='gray', lw=0.8, ls='--', alpha=0.5); ax.axvline(0.5, color='gray', lw=0.8, ls='--', alpha=0.5)
    ax.set_title(title, fontsize=11)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7.5, frameon=False)
    fig.text(0.01, -0.05, FOOTNOTE, fontsize=7.5, color='#555', ha='left')
    plt.tight_layout()
    save_fig(fig, out_stem)

def make_barplot(df, title, out_stem):
    scores = pd.to_numeric(df['Total'], errors='coerce').dropna().sort_values()
    fig, ax = plt.subplots(figsize=(8, max(3, len(scores) * 0.45)))
    bars = ax.barh(scores.index, scores.values, color=['orange' if 'BBKNN' in str(m) else 'steelblue' for m in scores.index])
    for bp, label in zip(bars, scores.index):
        if 'BBKNN' in str(label): bp.set_hatch('//')
        ax.text(bp.get_width() + 0.006, bp.get_y() + bp.get_height() / 2, f'{bp.get_width():.3f}', va='center', fontsize=8)
    ax.set_xlabel('Total aggregate score  (0.4 × Batch + 0.6 × Bio)', fontsize=10); ax.set_xlim(0, 1.08)
    ax.axvline(0.5, color='gray', ls='--', lw=0.8, alpha=0.5); ax.set_title(title, fontsize=11)
    fig.text(0.01, -0.05, FOOTNOTE, fontsize=7.5, color='#555', ha='left')
    plt.tight_layout()
    save_fig(fig, out_stem)

make_scatter(df_scaled, f'Batch correction vs Bio conservation (min-max)\n{PUB["name"]} · {TAG}', PUB['fig'])
make_results_table(df_scaled, f'scib — {PUB["name"]} · numbers min-max, color raw · {TAG}', PUB['tab'], color_df=df_fixed)
make_scatter(df_fixed, f'Batch correction vs Bio conservation (raw)\n{PUB["name"]} · {TAG}', 'scib_scatter_raw')
make_results_table(df_fixed, f'scib — {PUB["name"]} · raw · {TAG}', 'scib_results_table_raw')
make_barplot(df_scaled, f'scib Total (min-max) — {PUB["name"]} · {TAG}', 'scib_barplot_minmax')
print('wrote', OUT)
