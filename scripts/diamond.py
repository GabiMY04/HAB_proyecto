"""
Script: diamond.py

Descripción:
    Implementa el algoritmo DIAMOnD (Disease Module Detection) para expandir un conjunto
    de genes semilla dentro de una red génica. El método añade iterativamente los genes
    más significativamente conectados con las semillas mediante una prueba hipergeométrica.
    El algoritmo finaliza cuando se cumple alguno de los siguientes criterios:
        - Se alcanzan 200 iteraciones.
        - Se añaden 100 genes nuevos.
        - El valor p deja de ser significativo (p > 0.05).

Entradas:
    - network_*.txt           → red de entrada (lista de aristas).
    - genes_semilla_stringid.txt → lista de genes semilla (uno por línea).

Salidas:
    - diamond_results.csv → tabla con genes del módulo, p-valores y conexiones.
    - diamond_genes.txt   → lista de genes del módulo (uno por línea).
"""

from pathlib import Path
import pandas as pd
import networkx as nx
from scipy.stats import hypergeom


def leer_genes(input_genes: Path) -> list[str]:
    """Lee un archivo con genes y devuelve una lista limpia."""
    if not input_genes.exists():
        raise FileNotFoundError(f"No se encontró el archivo de entrada: {input_genes}")

    contenido = input_genes.read_text().strip()
    return [g.strip() for g in contenido.replace(",", "\n").splitlines() if g.strip()]


def ejecutar_diamond(input_dir: Path, output_dir: Path) -> pd.DataFrame:
    """Ejecuta el algoritmo DIAMOnD sobre una red génica."""
    # Configuración de rutas
    ruta_red = input_dir / "network_arabidopsis.txt"
    ruta_semillas = input_dir / "genes_semilla_stringid.txt"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Cargar red
    try:
        G = nx.read_edgelist(ruta_red, delimiter=None)
    except Exception:
        G = nx.read_edgelist(ruta_red, delimiter=" ")
    print(f"✅ Red cargada: {G.number_of_nodes()} nodos, {G.number_of_edges()} interacciones.")

    # Leer genes semilla
    semillas = set(leer_genes(ruta_semillas))

    # Normalizar IDs (eliminar sufijos tipo .1, .2)
    def limpiar_id(g):
        return g.split(".")[0] if g.count(".") > 1 else g

    semillas_limpias = {limpiar_id(g) for g in semillas}
    nodos_limpios = {limpiar_id(n) for n in G.nodes()}
    interseccion = semillas_limpias & nodos_limpios

    if not interseccion:
        print("⚠️ Ninguna semilla coincide directamente con la red. Intentando coincidencia parcial...")
        posibles = [n for n in G.nodes() if any(n.startswith(g) for g in semillas_limpias)]
        if not posibles:
            raise ValueError("Ninguna semilla se encuentra en la red, incluso tras normalización.")
        semillas = set(posibles)
    else:
        semillas = {n for n in G.nodes() if limpiar_id(n) in interseccion}

    print(f"✅ {len(semillas)} semillas válidas encontradas en la red.")
    print(f"Semillas iniciales ({len(semillas)}): {', '.join(sorted(semillas))}")

    K = len(semillas)
    N = G.number_of_nodes()

    # Precalcular vecinos y grados
    neighbors = {n: set(G.neighbors(n)) for n in G}
    degrees = {n: len(neighbors[n]) for n in G}

    def calc_pval(node):
        """Calcula p-valor hipergeométrico para un nodo."""
        k_i = degrees[node]
        k_s = len(neighbors[node] & semillas)
        return hypergeom.sf(k_s - 1, N, K, k_i), k_s

    # Inicializar candidatos
    candidatos = set().union(*[neighbors[g] for g in semillas]) - semillas
    results = []
    motivo_parada = None

    # Iteraciones principales
    for i in range(200):
        pvals = [(node, *calc_pval(node)) for node in candidatos]
        if not pvals:
            motivo_parada = "⚠️ Sin candidatos restantes"
            break

        node, pval, k_s = min(pvals, key=lambda x: x[1])
        if pval > 0.05:
            motivo_parada = f"🟡 P-valor no significativo (p={pval:.3e})"
            break

        semillas.add(node)
        results.append((node, pval, k_s))
        candidatos.update(neighbors[node])
        candidatos -= semillas

        if (i + 1) % 10 == 0 or i < 10:
            print(f"Iteración {i+1:03d}: añadido {node} (p={pval:.2e}, conexiones={k_s})")

        if len(results) >= 100:
            motivo_parada = "🔚 Alcanzado límite de 100 genes añadidos"
            break

    if motivo_parada is None:
        motivo_parada = "⏱️ Límite de 200 iteraciones alcanzado"

    # Guardar resultados
    df_added = pd.DataFrame(results, columns=["Gen", "p-valor", "Conexiones_con_semillas"])
    df_added["Tipo"] = "Añadido"

    semillas_restantes = semillas - set(df_added["Gen"])
    df_seeds = pd.DataFrame({
        "Gen": sorted(semillas_restantes),
        "p-valor": [None] * len(semillas_restantes),
        "Conexiones_con_semillas": [None] * len(semillas_restantes),
        "Tipo": ["Semilla"] * len(semillas_restantes),
    })

    df_final = pd.concat([df_seeds, df_added], ignore_index=True)
    df_final.attrs["criterio_parada"] = motivo_parada

    csv_path = output_dir / "diamond_results.csv"
    genes_path = output_dir / "diamond_genes.txt"
    df_final.to_csv(csv_path, index=False)

    # Escribir lista de genes
    genes_path.write_text("\n".join(df_final["Gen"].dropna()))

    print("\n✅ Algoritmo DIAMOnD completado.")
    print(f"Genes en módulo total: {len(df_final)} (Semillas + añadidos)")
    print(f"Criterio de parada: {motivo_parada}")
    print(f"Resultados guardados en:\n  • {csv_path}\n  • {genes_path}")

    return df_final


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Ejecuta el algoritmo DIAMOnD sobre una red génica.")
    parser.add_argument("--input_dir", type=Path, required=True, help="Directorio con los archivos de entrada (red y semillas).")
    parser.add_argument("--output_dir", type=Path, required=True, help="Directorio donde guardar los resultados.")
    args = parser.parse_args()

    ejecutar_diamond(args.input_dir, args.output_dir)