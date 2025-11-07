"""
Script: procesar_DEGs.py

Descripción:
    Filtra los genes diferencialmente expresados (DEGs) obtenidos de un experimento
    transcriptómico (por ejemplo, en Arabidopsis thaliana).
    Identifica genes significativamente sobreexpresados y subexpresados a partir
    de un archivo con columnas 'Gene', 'Coef' y 'P.value.adj'.

Entradas:
    - Archivo TSV con resultados de expresión diferencial.

Salidas:
    - genes_up.txt        → genes sobreexpresados.
    - genes_down.txt      → genes subexpresados.
    - genes_semilla.txt   → combinación de ambos conjuntos.
"""

import argparse
from pathlib import Path
import pandas as pd


def procesar_DEGs(input_file: Path, output_dir: Path, pval: float = 0.05) -> None:
    """Procesa los resultados de expresión diferencial y genera archivos con los genes significativos."""

    # Leer archivo TSV; si los IDs están en el índice, renombrarlo como 'Gene'
    df = pd.read_csv(input_file, sep="\t", index_col=0)
    df.index.name = "Gene"
    df.reset_index(inplace=True)

    # Verificar que existan las columnas necesarias
    columnas_requeridas = {"Gene", "Coef", "P.value.adj"}
    if not columnas_requeridas.issubset(df.columns):
        raise ValueError(f"El archivo debe contener las columnas: {', '.join(columnas_requeridas)}")

    # Filtrar genes significativos según p-valor ajustado
    df_sig = df[df["P.value.adj"] < pval]

    # Separar por dirección del cambio de expresión
    genes_up = df_sig.loc[df_sig["Coef"] > 0, "Gene"].tolist()
    genes_down = df_sig.loc[df_sig["Coef"] < 0, "Gene"].tolist()
    genes_semilla = genes_up + genes_down

    # Crear directorio de salida si no existe
    output_dir.mkdir(parents=True, exist_ok=True)

    # Guardar listas de genes
    (output_dir / "genes_up.txt").write_text("\n".join(genes_up))
    (output_dir / "genes_down.txt").write_text("\n".join(genes_down))
    (output_dir / "genes_semilla.txt").write_text("\n".join(genes_semilla))

    print(f"✅ {len(genes_up)} genes upregulados, {len(genes_down)} downregulados.")
    print(f"🧬 Total genes semilla: {len(genes_semilla)}")
    print(f"Archivos generados en: {output_dir}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Procesa los resultados de expresión diferencial y genera listas de genes.")
    parser.add_argument("--input_file", type=Path, required=True, help="Archivo TSV con resultados de expresión diferencial.")
    parser.add_argument("--output_dir", type=Path, default=Path("data/"), help="Directorio de salida.")
    parser.add_argument("--pval", type=float, default=0.05, help="Umbral de significación ajustada.")
    args = parser.parse_args()

    procesar_DEGs(args.input_file, args.output_dir, args.pval)