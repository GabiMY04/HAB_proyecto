"""
Script: visualizacion_omica.py

Descripción:
    Funciones utilitarias para generar visualizaciones ómicas a partir de:
    - Resultados de expresión diferencial (DEGs)
    - Redes de interacción (edgelist STRING / PPI)
    - Listas de genes semilla y genes añadidos por DIAMOnD

Funciones principales:
    - load_deg(path_pattern)
        Lee un archivo TSV de DEGs (acepta patrón glob) y devuelve un pandas.DataFrame.
        Formato recomendado del TSV: columnas mínimas 'Gene', 'Coef' (log2FC), 'P.value.adj'.

    - plot_volcano(df, coef_col='Coef', padj_col='P.value.adj', gene_col=None,
                   fc_thr=1.0, fdr_thr=0.05, out='results/plots/volcano.png', top_n=10)
        Genera un volcano plot (log2FC vs -log10(adj p-value)), marca genes significativos
        según umbrales y anota los top_n genes más significativos en la esquina superior derecha.
        Parámetros:
          * df: DataFrame con columnas de coeficiente y p-ajustado.
          * coef_col / padj_col: nombres de columnas en df.
          * out: ruta de salida PNG.
        Salida: imagen PNG guardada en la ruta indicada.

    - plot_network(edge_path, seeds_path=None, added_path=None,
                   out='results/plots/network_seed_overlay.png', show_labels=False)
        Dibuja un subgrafo relevante de la red PPI (semillas y añadidos resaltados).
        Entradas esperadas:
          * edge_path: edgelist (archivo plano, dos columnas: nodoA nodoB).
          * seeds_path / added_path: listas de genes (una línea por gen). Acepta IDs con o sin prefijo "3702.".
          * out: ruta de salida PNG.
        Salida: imagen PNG del grafo y (según versión) opcionalmente .graphml.


"""

import os
import glob
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy import stats
from networkx.algorithms.community import greedy_modularity_communities
import matplotlib.patches as mpatches


###################################################################
##################### load_deg  ###################################
###################################################################

def load_deg(path_pattern):
    """
    Carga el archivo de DEGs y devuelve un pandas.DataFrame
    """
    files = glob.glob(path_pattern)
    if not files:
        raise FileNotFoundError(f"No se encontró ningún archivo con patrón {path_pattern}")
    # toma el primero si hay varios
    df = pd.read_csv(files[0], sep=None, engine='python')
    return df



###################################################################
##################### plot_volcano ################################
###################################################################
def plot_volcano(df, coef_col='Coef', padj_col='P.value.adj', gene_col=None,
                 fc_thr=1.0, fdr_thr=0.05, out='results/plots/volcano.png', top_n=10):
    os.makedirs(os.path.dirname(str(out)), exist_ok=True)
    """
    Grafica un volcano plot y anota en la parte superior derecha los `top_n`
    genes más significativos (ordenados por mayor -log10(padj)).
    """
    df = df.copy()
    # Columnas de genes
    if gene_col is None:
        if 'Gene' in df.columns:
            gene_col = 'Gene'
        elif df.index.dtype == object:
            df = df.reset_index()
            gene_col = df.columns[0]
        else:
            gene_col = df.columns[0]

    df[padj_col] = pd.to_numeric(df[padj_col], errors='coerce')
    df[coef_col] = pd.to_numeric(df[coef_col], errors='coerce')
    df['neglog10padj'] = -np.log10(df[padj_col].replace(0, 1e-300))
    df['sig'] = (df[padj_col] < fdr_thr) & (df[coef_col].abs() >= fc_thr)

    plt.figure(figsize=(8,6))
    colors = df['sig'].map({True: 'red', False: 'gray'})
    plt.scatter(df[coef_col], df['neglog10padj'], c=colors, s=10, alpha=0.6)
    plt.axhline(-np.log10(fdr_thr), color='black', linestyle='--', lw=0.8)
    plt.axvline(-fc_thr, color='black', linestyle='--', lw=0.6)
    plt.axvline(fc_thr, color='black', linestyle='--', lw=0.6)
    plt.xlabel('log2FC (Coef)')
    plt.ylabel('-log10(adj p-value)')
    plt.title('Volcano plot')

    # preparar anotaciones
    top_genes = df.nlargest(top_n, 'neglog10padj').copy()
    if not top_genes.empty:
        ax = plt.gca()
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        x_range = x_max - x_min if (x_max - x_min) != 0 else 1.0
        y_range = y_max - y_min if (y_max - y_min) != 0 else 1.0

        label_x = x_max - 0.02 * x_range

        label_top = y_max - 0.02 * y_range
        label_bottom = y_max - 0.35 * y_range
        if label_bottom >= label_top:
            label_bottom = y_max - 0.15 * y_range


        title_y = min(y_max - 0.01 * y_range, label_top + 0.02 * y_range)
        ax.text(label_x, title_y, "Top 10", fontsize=9, fontweight='bold', ha='right', va='bottom')


        n = len(top_genes)
        label_ys = np.linspace(label_top - 0.02 * y_range, label_bottom, n)

        top_genes = top_genes.sort_values('neglog10padj', ascending=False).reset_index(drop=True)

        for i, row in top_genes.iterrows():
            px = float(row[coef_col])
            py = float(row['neglog10padj'])
            label = str(row[gene_col])
            y_label = label_ys[i]
            
            mid_x = px + 0.25 * (label_x - px)
            ax.plot([px, mid_x, label_x], [py, (py + y_label) / 2, y_label],
                    color='0.5', linewidth=0.7, solid_capstyle='round')

            ax.scatter([px], [py], c='black', s=20, zorder=5)

            ax.text(label_x, y_label, label, fontsize=8, ha='right', va='center',
                    bbox=dict(facecolor='white', alpha=0.0, edgecolor='none'))

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

    plt.tight_layout()
    plt.savefig(str(out), dpi=200)
    plt.close()
    print("Volcano guardado en", out)




###################################################################
##################### plot_network ################################
###################################################################
def plot_network(edge_path, seeds_path=None, added_path=None,
                 out='results/plots/network_seed_overlay.png',
                 show_labels=False):
    """
    Genera una visualización 
    de la red de interacciones proteína-proteína (PPI),
    destacando los genes semilla y los añadidos por DIAMOnD.
    """
    os.makedirs(os.path.dirname(out), exist_ok=True)

    # ===Cargar aristas ===
    edges = pd.read_csv(edge_path, sep=None, engine='python', header=None)
    if edges.shape[1] < 2:
        raise ValueError("Edge list debe tener al menos 2 columnas (nodeA nodeB)")
    G = nx.from_pandas_edgelist(edges, 0, 1)

    # === Listas de genes ===
    seeds = set()
    if seeds_path and os.path.exists(seeds_path):
        seeds = set(line.strip() for line in open(seeds_path) if line.strip())

    added = set()
    if added_path and os.path.exists(added_path):
        added = set(line.strip() for line in open(added_path) if line.strip())

    sub_nodes = set(seeds) | set(added)
    valid_nodes = set(G.nodes())

    for n in list(sub_nodes):
        if n in valid_nodes:
            sub_nodes |= set(G.neighbors(n))

    sub_nodes = sub_nodes & valid_nodes
    H = G.subgraph(sub_nodes).copy()

    # === Detectar comunidades ===
    try:
        communities = list(greedy_modularity_communities(H))
        color_map = {}
        for i, com in enumerate(communities):
            for node in com:
                color_map[node] = i
        community_colors = [color_map[n] for n in H.nodes()]
    except Exception:
        community_colors = ["lightgray"] * len(H.nodes())


    pos = nx.spring_layout(H, seed=42, k=0.45, iterations=150)
    plt.figure(figsize=(13, 11))
    plt.axis('off')

    # === Aristas ===
    nx.draw_networkx_edges(H, pos, alpha=0.15, width=0.5, edge_color="gray")

    # === Nodos de fondo ===
    nx.draw_networkx_nodes(
        H, pos,
        node_color='lightgray',
        node_size=[20 + 4 * H.degree(n) for n in H.nodes()],
        alpha=0.5,
        edgecolors='none'
    )

    # ===Añadidos (naranja) ===
    added_nodes = list(added & set(H.nodes()))
    if added_nodes:
        nx.draw_networkx_nodes(
            H, pos,
            nodelist=added_nodes,
            node_color='orange',
            node_size=150,
            edgecolors='k',
            linewidths=0.4
        )

    # === Semillas (rojo) ===
    seed_nodes = list(seeds & set(H.nodes()))
    if seed_nodes:
        nx.draw_networkx_nodes(
            H, pos,
            nodelist=seed_nodes,
            node_color='red',
            node_size=220,
            edgecolors='k',
            linewidths=0.4
        )

    if show_labels:
        labels = {n: n for n in H.nodes() if n in seeds or n in added}
        nx.draw_networkx_labels(H, pos, labels=labels, font_size=6)

    # === Leyenda
    red_patch = mpatches.Circle((0, 0), radius=0.1, color='red', label='Genes semilla')
    orange_patch = mpatches.Circle((0, 0), radius=0.1, color='orange', label='Genes añadidos (DIAMOnD)')
    gray_patch = mpatches.Circle((0, 0), radius=0.1, color='lightgray', label='Otros genes')

    legend_elements = [red_patch, orange_patch, gray_patch]
    plt.legend(
        handles=legend_elements,
        loc='upper right',
        fontsize=10,
        frameon=True,
        facecolor='white',
        framealpha=0.9,
        borderpad=0.8,
        title="Leyenda",
        title_fontsize=11
    )

    plt.title(
        "Red de Interacciones Proteína–Proteína (PPI)\n",
        fontsize=13,
        pad=15
    )

    # === Guardar figura ===
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()

    graphml_path = out.replace(".png", ".graphml")
    nx.write_graphml(H, graphml_path)
