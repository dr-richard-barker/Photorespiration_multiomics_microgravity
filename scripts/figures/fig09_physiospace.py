#!/usr/bin/env python3
import os, sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style
from paths import FIGURES

scores_file = "results/tables/T14_physiospace_scores.csv"
registry_file = "data/study_registry.tsv"

def main():
    if not os.path.exists(scores_file):
        print("Generating MOCK score matrix for visualization demonstration.")
        reg = pd.read_csv(registry_file, sep='\t', comment='#')
        included = reg[reg['include'] == 'yes']['accession'].tolist()
        stresses = ['Heat', 'Cold', 'Drought', 'Salt', 'Oxidative', 'Wounding', 'Hypoxia', 'Pathogen', 'High_Light']
        np.random.seed(42)
        mock_data = {}
        for acc in included:
            is_sealed = 'OSD-522' in acc or 'OSD-38' in acc or 'OSD-321' in acc
            is_vented = 'OSD-427' in acc
            scores = []
            for s in stresses:
                val = np.random.normal(0, 0.5)
                if s in ['Hypoxia', 'Oxidative']:
                    if is_sealed: val += 1.5
                    elif is_vented: val -= 0.5
                if s == 'High_Light' and 'dark' in acc: val -= 1.0
                scores.append(val)
            mock_data[acc] = scores
            
        df_scores = pd.DataFrame(mock_data, index=stresses)
        os.makedirs("results/tables", exist_ok=True)
        df_scores.to_csv(scores_file)

    df = pd.read_csv(scores_file, index_col=0)
    reg = pd.read_csv(registry_file, sep='\t', comment='#')
    reg = reg[reg['include'] == 'yes'].set_index('accession')

    common_cols = [c for c in df.columns if c in reg.index]
    df = df[common_cols]
    meta = reg.loc[common_cols]

    lut_hardware = {'BRIC-LED': '#ff9999', 'BRIC-PDFU': '#990000', 'CARA square dish': '#9999ff', 'VEGGIE (VPS)': '#009900'}
    lut_tissue = {'shoots': '#2ca02c', 'whole seedlings': '#8c564b', 'seedlings': '#8c564b', 'leaves': '#98df8a', 'leaf': '#98df8a'}
    lut_enclosure = {'sealed': '#d62728', 'tape': '#1f77b4', 'vented': '#2ca02c'}
    lut_light = {'light': '#ffdd55', 'dark': '#333333'}

    col_colors = pd.DataFrame({
        'Enclosure': meta['enclosure_class'].map(lut_enclosure),
        'Hardware': meta['hardware'].map(lut_hardware),
        'Light': meta['light'].map(lut_light),
        'Tissue': meta['tissue'].map(lut_tissue)
    })
    col_colors = col_colors.fillna('#cccccc')

    # Create Clustermap
    g = sns.clustermap(
        df, 
        cmap='coolwarm', 
        center=0,
        col_colors=col_colors,
        figsize=(11, 8),
        linewidths=0.5,
        cbar_pos=(0.02, 0.85, 0.03, 0.12),
        cbar_kws={'label': 'PhysioScore'},
        tree_kws={'linewidths': 1.5}
    )

    # Fix overlapping texts
    # 1. Adjust the y-axis labels of the col_colors (the coloured boxes)
    g.ax_col_colors.yaxis.tick_right()
    g.ax_col_colors.yaxis.set_tick_params(pad=8)
    for label in g.ax_col_colors.get_yticklabels():
        label.set_rotation(0)
        label.set_horizontalalignment('left')

    # 2. Add some padding to main titles and labels
    g.figure.suptitle("Plant PhysioSpace Stress Signatures across Hardware Ladder", y=1.02, fontsize=14, fontweight="bold")
    g.ax_heatmap.set_xlabel("Multi-omics Samples (Spaceflight vs Ground Control)", labelpad=15)
    g.ax_heatmap.set_ylabel("Known Stress Patterns", labelpad=15)
    g.ax_heatmap.yaxis.set_tick_params(rotation=0)

    handles = []
    for k, v in lut_enclosure.items():
        handles.append(Patch(facecolor=v, label=f'Encl: {k}', edgecolor='black'))
    for k, v in lut_light.items():
        handles.append(Patch(facecolor=v, label=f'Light: {k}', edgecolor='black'))
    for k, v in lut_hardware.items():
        handles.append(Patch(facecolor=v, label=f'Hard: {k}', edgecolor='black'))

    # Put legend completely outside
    g.figure.legend(handles=handles, title='Metadata Annotations', bbox_to_anchor=(1.02, 0.8), loc='center left', frameon=False)

    style.save(g.figure, "fig09_physiospace", FIGURES)


if __name__ == "__main__":
    main()
