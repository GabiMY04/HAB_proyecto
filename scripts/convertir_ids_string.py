"""
Script: convertir_ids_string.py

Descripción:
    Convierte una lista de genes de Arabidopsis thaliana en formato AGI
    (por ejemplo, AT1G54410) al formato utilizado por STRINGdb
    (por ejemplo, 3702.AT1G54410), agregando el prefijo taxonómico '3702'.

Entradas:
    - Lista de genes en formato AGI (uno por línea).

Salidas:
    - Lista de genes con prefijo STRINGdb (3702.<GENE>).
"""

import argparse
from pathlib import Path


def convertir_a_string_ids(input_file: Path, output_file: Path) -> None:
    """Convierte los identificadores AGI a formato STRINGdb (3702.<GENE>)."""

    # Leer genes válidos del archivo
    genes = [line.strip() for line in input_file.read_text().splitlines() if line.strip()]

    # Validar formato básico AGI (ATxGxxxxx)
    genes_validos = [g for g in genes if g.startswith("AT") and "G" in g]
    if len(genes_validos) < len(genes):
        print(f"⚠️ {len(genes) - len(genes_validos)} genes no válidos fueron omitidos.")

    # Agregar prefijo taxonómico de Arabidopsis (3702)
    genes_string = [f"3702.{g}" for g in genes_validos]

    # Crear directorio de salida si no existe
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Guardar resultado
    output_file.write_text("\n".join(genes_string))
    print(f"✅ {len(genes_string)} genes convertidos al formato STRINGdb.")
    print(f"📁 Archivo guardado en: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convierte identificadores AGI de Arabidopsis al formato STRINGdb.")
    parser.add_argument("--input_file", type=Path, required=True, help="Archivo con genes en formato AGI (uno por línea).")
    parser.add_argument("--output_file", type=Path, default=Path("data/genes_semilla_string.txt"), help="Archivo de salida.")
    args = parser.parse_args()

    convertir_a_string_ids(args.input_file, args.output_file)