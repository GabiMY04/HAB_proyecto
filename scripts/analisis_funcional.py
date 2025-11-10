"""
Script: analisis_funcional.py

Descripción:
    Realiza un análisis de sobrerrepresentación génica (Over-Representation Analysis, ORA)
    utilizando la API oficial de STRINGdb. Permite analizar listas de genes de Arabidopsis
    thaliana y devuelve solo los términos de Gene Ontology (GO), reemplazando los
    identificadores GO:XXXXXXX por sus nombres descriptivos.

Entradas:
    - Archivo con genes (uno por línea o separados por comas).

Salidas:
    - enrichment_results.csv → categorías GO con nombre descriptivo y FDR.
    - enrichment_plot.png    → gráfico de barras con las categorías más significativas.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import requests
from goatools.obo_parser import GODag


# === Utilidades ===
def leer_genes(input_genes: str) -> list[str]:
    """Lee un archivo con genes (separados por comas o saltos de línea) y devuelve una lista."""
    path = Path(input_genes)
    if not path.exists():
        raise FileNotFoundError(f"No se encontró el archivo de entrada: {path}")

    contenido = path.read_text().strip()
    genes = [g.strip() for g in contenido.replace(",", "\n").splitlines() if g.strip()]
    print(f"✅ {len(genes)} IDs leídos desde {path}")
    return genes


def mapear_GO_a_nombre(go_ids: list[str], output_dir: Path = Path(".")) -> dict[str, str]:
    """Mapea identificadores GO (GO:XXXXXXX) a nombres descriptivos usando la ontología GO."""
    obo_path = output_dir / "go-basic.obo"

    # Descargar la ontología si no existe
    if not obo_path.exists():
        print("⬇️ Descargando ontología GO básica (go-basic.obo)...")
        url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
        r = requests.get(url)
        r.raise_for_status()
        obo_path.write_bytes(r.content)

    go_dag = GODag(obo_path, optional_attrs={"defn"})
    return {go_id: go_dag[go_id].name for go_id in go_ids if go_id in go_dag}

def obtener_enriquecimiento_STRING(genes: list[str], species: int = 3702) -> pd.DataFrame:
    """
    Consulta la API oficial de STRINGdb para realizar un análisis de sobrerrepresentación (ORA).
    Devuelve un DataFrame con los resultados de enriquecimiento.
    """
    url = "https://string-db.org/api/json/enrichment"
    params = {"identifiers": "%0d".join(genes), "species": species}
    r = requests.get(url, params=params)
    r.raise_for_status()
    data = r.json()
    if not data:
        return pd.DataFrame()
    return pd.DataFrame(data)


def graficar_resultados(df: pd.DataFrame, output_dir: Path, n_resultados: int = 10) -> None:
    """
    Genera un gráfico de barras con los n_resultados más significativos del análisis ORA.

    Los valores de significancia se representan como -log10(FDR), lo que permite visualizar
    de forma intuitiva la fuerza del enriquecimiento. Se emplean únicamente categorías GO
    ya mapeadas a descripciones legibles.
    """
    # Detectar la columna de p-valor ajustado
    if "FDR" in df.columns:
        pvals = "FDR"
        label = "-log10(FDR)"
    elif "fdr" in df.columns:
        pvals = "fdr"
        label = "-log10(FDR)"
    elif "Adjusted P-value" in df.columns:
        pvals = "Adjusted P-value"
        label = "-log10(Adjusted P-value)"
    else:
        raise KeyError("No se encontró una columna válida de p-valores ajustados (FDR).")

    # Ordenar por significancia y seleccionar top N
    df = df.sort_values(pvals).head(n_resultados).iloc[::-1].copy()

    # Acortar nombres largos
    df["Category_short"] = df["Category"].apply(lambda t: t[:67] + "..." if isinstance(t, str) and len(t) > 70 else t)

    # Crear gráfico
    plt.figure(figsize=(10, 0.6 * len(df)))
    plt.barh(df["Category_short"], -np.log10(df[pvals]), color="#1976D2")
    plt.xlabel(label, fontsize=12)
    plt.title("Categorías GO más representadas (STRINGdb)", fontsize=14, pad=15)
    plt.tight_layout()

    # Guardar gráfico
    output_png = output_dir / "enrichment_plot.png"
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close()


# === ORA principal ===
def ejecutar_ora_STRING(input_genes: str, output_dir: str) -> pd.DataFrame:
    """
    Ejecuta un análisis ORA en STRINGdb y devuelve los resultados limitados a términos GO.
    Los identificadores GO se reemplazan por sus descripciones oficiales.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Leer genes
    genes = leer_genes(input_genes)
    print(f"🔍 Ejecutando ORA en STRINGdb para {len(genes)} genes (species={3702})...")

    # Ejecutar análisis ORA
    try:
        enrichment = obtener_enriquecimiento_STRING(genes, species=3702)
        if enrichment is None or enrichment.empty:
            print("⚠️ No se obtuvieron resultados de enriquecimiento.")
            return pd.DataFrame()
    except Exception as e:
        raise RuntimeError(f"❌ Error al ejecutar enriquecimiento en STRINGdb: {e}")

    # Filtrar solo términos GO
    mask_go = enrichment["term"].str.startswith(("GO:", "GOBP:", "GOCC:", "GOMF:"))
    enrichment = enrichment[mask_go].copy()
    if enrichment.empty:
        print("⚠️ No se encontraron términos GO en los resultados.")
        return pd.DataFrame()

    # Mapeo GO → nombre descriptivo
    print("ℹ️ Mapeando identificadores GO a nombres descriptivos...")
    go_ids = enrichment["term"].unique().tolist()
    mapping = mapear_GO_a_nombre(go_ids, output_dir)
    enrichment["Category"] = enrichment["term"].map(mapping).fillna(enrichment["description"])

    # Seleccionar columnas relevantes
    columnas_exportar = [
        "Category",
        "category",
        "fdr" if "fdr" in enrichment.columns else "p_value",
        "inputGenes" if "inputGenes" in enrichment.columns else "number_of_genes",
    ]
    columnas_exportar = [c for c in columnas_exportar if c in enrichment.columns]
    df_final = enrichment[columnas_exportar].rename(columns={"fdr": "FDR"})

    # Guardar resultados
    output_csv = output_dir / "enrichment_results.csv"
    df_final.to_csv(output_csv, index=False)
    print(f"✅ Resultados guardados en {output_csv.parent}/")

    # Graficar resultados
    try:
        graficar_resultados(df_final, output_dir, 10)
    except Exception:
        print("⚠️ No se pudo generar el gráfico, pero los resultados fueron guardados correctamente.")

    return df_final




# =============== Propagación GUILD (Random Walk with Restart) =====================


from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
from scipy import sparse

def _leer_lista(path: Path) -> list[str]:
    return [x.strip() for x in Path(path).read_text().splitlines() if x.strip()]

def ejecutar_guild(input_dir: Path, out_propagation: Path, topk: int = 100, alpha: float = 0.5,
                   tol: float = 1e-8, max_iter: int = 200) -> pd.DataFrame:
    """
    Ejecuta GUILD (priorización por propagación en red; RWR) sobre la red y semillas.
    - input_dir debe contener: network_arabidopsis.txt y genes_semilla_stringid.txt
    - out_propagation: carpeta donde guardar guild_scores.csv y guild_genes.txt
    """
    out_propagation.mkdir(parents=True, exist_ok=True)
    ruta_red = input_dir / "network_arabidopsis.txt"
    ruta_sem = input_dir / "genes_semilla_stringid.txt"

    # Cargar red (igual patrón que DIAMOnD)
    try:
        G = nx.read_edgelist(ruta_red, delimiter=None)
    except Exception:
        G = nx.read_edgelist(ruta_red, delimiter=" ")

    semillas = set(_leer_lista(ruta_sem))
    if not semillas:
        raise ValueError("No hay semillas en genes_semilla_stringid.txt")

    # Normalizar IDs como en DIAMOnD (quitar sufijos .1, .2...)
    def _clean(s): return s.split(".")[0] if s.count(".") > 1 else s
    semillas_norm = {_clean(s) for s in semillas}
    nodos_norm = {_clean(n) for n in G.nodes()}
    inter = semillas_norm & nodos_norm
    if inter:
        semillas = {n for n in G.nodes() if _clean(n) in inter}
    else:
        posibles = [n for n in G.nodes() if any(n.startswith(s) for s in semillas_norm)]
        if not posibles:
            raise ValueError("Ninguna semilla está en la red tras normalización.")
        semillas = set(posibles)

    # Construir matriz columna-normalizada W
    nodes = list(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    rows, cols = [], []
    for u, v in G.edges():
        rows += [idx[v]]; cols += [idx[u]]
        rows += [idx[u]]; cols += [idx[v]]
    data = np.ones(len(rows), dtype=float)
    A = sparse.coo_matrix((data, (rows, cols)), shape=(len(nodes), len(nodes))).tocsr()
    col_sums = np.array(A.sum(axis=0)).ravel()
    col_sums[col_sums == 0] = 1.0
    W = A @ sparse.diags(1.0 / col_sums)

    # Vector de reinicio r
    r = np.zeros(len(nodes))
    seeds_idx = [idx[s] for s in semillas if s in idx]
    r[seeds_idx] = 1.0 / len(seeds_idx)

    # RWR: p_{t+1} = (1 - alpha) * W * p_t + alpha * r
    p = r.copy()
    for _ in range(max_iter):
        p_new = (1 - alpha) * (W @ p) + alpha * r
        if np.linalg.norm(p_new - p, 1) < tol:
            p = p_new; break
        p = p_new

    # Tabla de scores y selección top-K no semillas
    df_scores = pd.DataFrame({"Gen": nodes, "Score": p})
    df_scores = df_scores.sort_values("Score", ascending=False)
    df_scores.to_csv(out_propagation / "guild_scores.csv", index=False)

    candidatos = df_scores[~df_scores["Gen"].isin(semillas)].head(topk)["Gen"].tolist()
    genes_modulo = list(semillas) + candidatos
    (out_propagation / "guild_genes.txt").write_text("\n".join(genes_modulo))

    print(f"✅ GUILD completado. {len(candidatos)} genes añadidos.")
    print(f"   • {out_propagation/'guild_scores.csv'}")
    print(f"   • {out_propagation/'guild_genes.txt'}")
    return df_scores
