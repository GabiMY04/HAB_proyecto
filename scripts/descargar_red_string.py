"""
Script: descargar_red_string.py

Descripción:
    Descarga una red de interacciones proteína-proteína (PPI) desde STRINGdb
    centrada en un conjunto de genes semilla. Filtra las interacciones según
    un umbral mínimo de score y guarda la red en formato de lista de aristas,
    compatible con NetworkX y DIAMOnD.

Entradas:
    - Archivo de texto con identificadores de genes (uno por línea).

Salidas:
    - genes_semilla_stringid.txt   → genes mapeados a STRING IDs.
    - network_subgraph.txt         → red de interacciones filtrada.
"""

import argparse
from pathlib import Path
import requests
import pandas as pd


def mapear_genes_a_STRING(genes: list[str], species_id: int) -> list[str]:
    """Convierte una lista de genes al formato de STRING usando su API oficial."""
    url = "https://string-db.org/api/json/get_string_ids"
    ids_mapeados = []

    for gene in genes:
        base = gene.replace(f"{species_id}.", "")  # quitar prefijo previo
        r = requests.get(url, params={"identifiers": base, "species": species_id}, timeout=10)
        if r.status_code == 200:
            data = r.json()
            if data:
                ids_mapeados.append(data[0]["stringId"])
                continue
        # si no hay coincidencia exacta
        ids_mapeados.append(f"{species_id}.{base}.1")

    return ids_mapeados


def descargar_red_STRING(genes_file: Path, species_id: int, output_file: Path, score_threshold: int = 700) -> None:
    """
    Descarga una red de interacciones proteína-proteína (PPI) desde la API de STRINGdb,
    centrada en los genes semilla especificados.

    El proceso garantiza compatibilidad entre los identificadores de los genes y los nodos
    de la red, filtrando interacciones por score mínimo. Produce una lista de aristas
    lista para su uso en NetworkX o DIAMOnD.
    """

    # Leer genes desde el archivo
    genes = [g.strip() for g in genes_file.read_text().splitlines() if g.strip()]
    print(f"🧬 Mapeando {len(genes)} genes a STRING IDs válidos...")

    # Mapear a IDs de STRING
    string_ids = mapear_genes_a_STRING(genes, species_id)

    # Guardar lista mapeada
    map_file = genes_file.parent / "genes_semilla_stringid.txt"
    map_file.write_text("\n".join(string_ids))
    print(f"📁 Mapeo guardado en: {map_file}")

    # Descargar red centrada en los genes
    identifiers = "%0d".join(string_ids)
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": identifiers,
        "species": species_id,
        "required_score": score_threshold,
        "add_nodes": 500,
        "network_type": "functional",
    }

    print(f"⏳ Descargando red centrada en las semillas (species={species_id}, score≥{score_threshold})...")
    r = requests.get(url, params=params)
    if r.status_code != 200:
        raise ConnectionError(f"Error {r.status_code}: {r.text}")

    data = r.json()
    if not data:
        raise ValueError("STRING no devolvió resultados. Verifique los genes o el species_id.")

    # Guardar lista de aristas (sin encabezado)
    df = pd.DataFrame(data)
    df[["stringId_A", "stringId_B"]].to_csv(output_file, sep=" ", index=False, header=False)

    print(f"✅ Red guardada en: {output_file}")
    print(f"🔗 Total de interacciones: {len(df)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Descarga una red de interacciones centrada en genes semilla desde STRINGdb.")
    parser.add_argument("--genes_file", type=Path, required=True, help="Archivo con genes semilla (uno por línea).")
    parser.add_argument("--species_id", type=int, required=True, help="Código taxonómico de la especie (ej. 3702, 9606, etc.).")
    parser.add_argument("--output_file", type=Path, default=Path("data/network_subgraph.txt"), help="Archivo de salida.")
    parser.add_argument("--score_threshold", type=int, default=700, help="Umbral mínimo de score combinado (0–1000).")
    args = parser.parse_args()

    descargar_red_STRING(args.genes_file, args.species_id, args.output_file, args.score_threshold)