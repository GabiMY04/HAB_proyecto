"""
Script: generar_visualizaciones.py

Descripción:
    Genera comparaciones gráficas de los resultados del análisis de enriquecimiento funcional
    antes y después de aplicar el algoritmo DIAMOnD. Funciona con resultados obtenidos desde
    STRINGdb y produce:
        - Comparación general de significancia (Pre vs Post)
        - Categorías nuevas tras la propagación
        - Diagrama de solapamiento (Venn)

Entradas:
    - results/ORA_semillas/enrichment_results.csv   → resultados pre-DIAMOnD.
    - results/ORA_diamond/enrichment_results.csv    → resultados post-DIAMOnD.

Salidas:
    - ora_barplot_delta.png   → comparación general Pre vs Post.
    - ora_barplot_nuevas.png  → categorías nuevas tras DIAMOnD.
    - ora_venn.png            → solapamiento entre categorías.
"""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


def generar_visualizaciones(pre_csv: Path, post_csv: Path, output_dir: Path) -> pd.DataFrame:
    """
    Genera representaciones visuales comparando los resultados del análisis funcional
    antes y después de la propagación con DIAMOnD.

    Incluye:
        - Comparación de significancia (barplot delta)
        - Nuevas categorías enriquecidas post-propagación
        - Diagrama de Venn del solapamiento de términos GO
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pre = pd.read_csv(pre_csv)
    post = pd.read_csv(post_csv)

    # Normalizar nombres de columnas
    for df in (pre, post):
        df.columns = df.columns.str.strip().str.lower().str.replace(" ", "_")
        df = df.loc[:, ~df.columns.duplicated()]

    # Mantener solo 'category' si coexiste con 'term'
    if "category" in pre.columns and "term" in pre.columns:
        pre = pre.drop(columns=["term"], errors="ignore")
    if "category" in post.columns and "term" in post.columns:
        post = post.drop(columns=["term"], errors="ignore")

    # Determinar columna principal de término
    if "category" in pre.columns:
        term_col = "category"
    elif "term" in pre.columns:
        term_col = "term"
    else:
        raise KeyError("❌ No se encontró ninguna columna válida ('category' o 'term') en los archivos ORA.")

    # Detectar columnas de p-valor
    p_pre_col = next((c for c in pre.columns if c in ("fdr", "adjusted_p-value", "adj_p_value", "p_value")), None)
    p_post_col = next((c for c in post.columns if c in ("fdr", "adjusted_p-value", "adj_p_value", "p_value")), None)
    if not p_pre_col or not p_post_col:
        raise KeyError("❌ No se encontraron columnas de p-valor ajustado ('fdr', 'adjusted_p-value', etc.).")

    # Seleccionar columnas relevantes
    pre = pre[[term_col, p_pre_col]].rename(columns={term_col: "term", p_pre_col: "p_pre"})
    post = post[[term_col, p_post_col]].rename(columns={term_col: "term", p_post_col: "p_post"})

    # Eliminar duplicados antes del merge
    pre = pre.loc[:, ~pre.columns.duplicated()]
    post = post.loc[:, ~post.columns.duplicated()]

    # Combinar datasets
    combined = pd.merge(pre, post, on="term", how="outer", suffixes=("_pre", "_post"))
    combined["estado"] = combined.apply(
        lambda r: (
            "compartido" if pd.notna(r["p_pre"]) and pd.notna(r["p_post"])
            else "solo_pre" if pd.notna(r["p_pre"])
            else "solo_post"
        ),
        axis=1,
    )

    combined["logp_pre"] = -np.log10(combined["p_pre"].fillna(1))
    combined["logp_post"] = -np.log10(combined["p_post"].fillna(1))
    combined["delta_logp"] = combined["logp_post"] - combined["logp_pre"]

    # === 1. Barplot comparativo general ===
    top_shared = combined[combined["estado"] == "compartido"].nsmallest(10, "p_post")
    plt.figure(figsize=(10, 0.6 * len(top_shared)))

    y = np.arange(len(top_shared))
    width = 0.35
    plt.barh(y - width / 2, top_shared["logp_pre"], width, label="Pre-DIAMOnD", color="#64B5F6")
    plt.barh(y + width / 2, top_shared["logp_post"], width, label="Post-DIAMOnD", color="#1976D2")

    plt.yticks(y, top_shared["term"], fontsize=9)
    plt.xlabel("-log10(FDR)")
    plt.title("Comparación del enriquecimiento funcional (Pre vs Post DIAMOnD)", fontsize=13, pad=15)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "ora_barplot_delta.png", dpi=300)
    plt.close()

    # === 2. Barplot de categorías nuevas ===
    new = combined[combined["estado"] == "solo_post"].nsmallest(10, "p_post")

    def acortar_texto(texto: str, max_len: int = 70) -> str:
        if isinstance(texto, str) and len(texto) > max_len:
            return texto[:max_len - 3] + "..."
        return texto

    new["term_short"] = new["term"].apply(acortar_texto)

    plt.figure(figsize=(8, 5 + 0.4 * len(new)))
    plt.barh(new["term_short"], -np.log10(new["p_post"]), color="#1E88E5")
    plt.title("Categorías nuevas tras DIAMOnD", fontsize=13)
    plt.xlabel("-log10(FDR)")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(output_dir / "ora_barplot_nuevas.png", dpi=300, bbox_inches="tight")
    plt.close()

    # === 3. Diagrama de Venn ===
    pre_terms = set(pre["term"])
    post_terms = set(post["term"])
    plt.figure(figsize=(5, 5))
    venn2([pre_terms, post_terms], set_labels=("Pre-DIAMOnD", "Post-DIAMOnD"))
    plt.title("Solapamiento de categorías significativas", fontsize=12)
    plt.tight_layout()
    plt.savefig(output_dir / "ora_venn.png", dpi=300)
    plt.close()

    print("✅ Gráficos generados correctamente:")
    print(f"   • {output_dir}/ora_barplot_delta.png")
    print(f"   • {output_dir}/ora_barplot_nuevas.png")
    print(f"   • {output_dir}/ora_venn.png")

    return combined


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Genera comparaciones gráficas de ORA pre y post DIAMOnD (STRINGdb).")
    parser.add_argument("--pre", type=Path, required=True, help="Archivo CSV del ORA pre-DIAMOnD.")
    parser.add_argument("--post", type=Path, required=True, help="Archivo CSV del ORA post-DIAMOnD.")
    parser.add_argument("--output_dir", type=Path, required=True, help="Directorio de salida para los gráficos.")

    args = parser.parse_args()
    generar_visualizaciones(args.pre, args.post, args.output_dir)