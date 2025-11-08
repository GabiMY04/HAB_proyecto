"""
Script: ejecutar_pipeline.py

Descripción:
    Orquesta el flujo completo de análisis funcional y propagación en red a partir
    de resultados de expresión diferencial (DEGs). Integra los pasos de procesamiento
    de DEGs, conversión de identificadores, descarga o reutilización de la red STRING,
    análisis funcional pre y post propagación (ORA), ejecución del algoritmo DIAMOnD,
    y generación de gráficas comparativas opcionales.

Entradas:
    - Archivo de expresión diferencial (TSV) con columnas: Gene, Coef, P.value.adj.
    - Archivos intermedios generados automáticamente:
        • data/genes_semilla.txt
        • data/genes_semilla_string.txt
        • data/genes_semilla_stringid.txt
        • data/network_arabidopsis.txt

Salidas:
    - results/ORA_semillas/enrichment_results.csv
    - results/diamond_propagation/diamond_results.csv
    - results/ORA_diamond/enrichment_results.csv
    - results/comparativas/*.png (opcional)
"""

import argparse
import warnings
import logging
from pathlib import Path

from scripts.procesar_DEGs import procesar_DEGs
from scripts.descargar_red_string import descargar_red_STRING
from scripts.convertir_ids_string import convertir_a_string_ids
from scripts.analisis_funcional import ejecutar_ora_STRING
from scripts.diamond import ejecutar_diamond
from scripts.comparar_enriquecimientos import generar_visualizaciones

# Silenciar advertencias y logs no relevantes
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
logging.getLogger("goatools").setLevel(logging.ERROR)

# Archivo de entrada principal
DEGS_RESULT_FILENAME = "Allcontrasts_GLM-Treat_P-0.1_FC-1.25_2025-10-14_16.57.27.tsv"


def main():
    parser = argparse.ArgumentParser(
        description="Pipeline de análisis funcional y propagación en red (DIAMOnD + STRINGdb)."
    )
    parser.add_argument("--input_dir", default="data/", help="Directorio con los archivos de entrada.")
    parser.add_argument("--output_dir", default="results/", help="Directorio donde se guardarán los resultados.")
    parser.add_argument("--species_id", type=int, default=3702, help="ID de especie en STRING (ej. 3702 = Arabidopsis thaliana).")
    args = parser.parse_args()

    data_dir = Path(args.input_dir)
    results_dir = Path(args.output_dir)
    DEGs_result_path = data_dir / DEGS_RESULT_FILENAME

    if not DEGs_result_path.exists():
        raise FileNotFoundError(f"No se encontró el archivo de resultados: {DEGs_result_path}")

    paso = 1
    print("🚀 Iniciando flujo de trabajo completo")

    # === Paso 1: Procesamiento de DEGs ===
    print(f"\n🔹 Paso {paso}: Procesando DEGs y generando lista de genes semilla...")
    procesar_DEGs(input_file=DEGs_result_path, output_dir=data_dir)
    paso += 1

    genes_semilla = data_dir / "genes_semilla.txt"
    if not genes_semilla.exists():
        raise FileNotFoundError(f"No se generó el archivo esperado: {genes_semilla}")

    # === Paso 2: Conversión a IDs de STRING ===
    print(f"\n🔹 Paso {paso}: Convirtiendo genes a formato STRING...")
    genes_semilla_string = data_dir / "genes_semilla_string.txt"
    convertir_a_string_ids(input_file=genes_semilla, output_file=genes_semilla_string)
    paso += 1

    # === Paso 3: Descarga o reutilización de red ===
    print(f"\n🔹 Paso {paso}: Preparando red de interacción...")
    input_network = data_dir / "network_arabidopsis.txt"

    if input_network.exists():
        print(f"🟢 Red existente detectada: {input_network.name}. Se reutilizará.")
    else:
        print("⏳ Red no encontrada. Descargando desde STRING...")
        descargar_red_STRING(
            genes_file=genes_semilla_string,
            output_file=input_network
        )
        print("✅ Red descargada correctamente.")
    paso += 1

    # === Paso 4: ORA pre-propagación ===
    print(f"\n🔹 Paso {paso}: Análisis funcional (ORA) de los genes semilla...")
    ora_pre_dir = results_dir / "ORA_semillas"
    genes_semilla_string_id = data_dir / "genes_semilla_stringid.txt"
    ejecutar_ora_STRING(
        input_genes=genes_semilla_string_id,
        output_dir=ora_pre_dir
    )
    paso += 1

    # === Paso 5: Propagación con DIAMOnD ===
    print(f"\n🔹 Paso {paso}: Propagación en red mediante DIAMOnD...")
    diamond_results_dir = results_dir / "diamond_propagation"
    ejecutar_diamond(input_dir=data_dir, output_dir=diamond_results_dir)
    paso += 1

    diamond_genes_path = diamond_results_dir / "diamond_genes.txt"
    if not diamond_genes_path.exists():
        raise FileNotFoundError(f"No se generó el archivo esperado: {diamond_genes_path}")

    # === Paso 6: ORA post-propagación ===
    print(f"\n🔹 Paso {paso}: Análisis funcional (ORA) de los genes propagados...")
    ora_post_dir = results_dir / "ORA_diamond"
    ejecutar_ora_STRING(
        input_genes=diamond_genes_path,
        output_dir=ora_post_dir,
        species_id=args.species_id
    )
    paso += 1

    # === Verificar resultados ===
    ora_pre_result = ora_pre_dir / "enrichment_results.csv"
    ora_post_result = ora_post_dir / "enrichment_results.csv"

    if not ora_pre_result.exists() or not ora_post_result.exists():
        raise FileNotFoundError(
            "❌ No se encontraron los resultados del análisis funcional (ORA).\n"
            f"Esperados:\n - {ora_pre_result}\n - {ora_post_result}"
        )

    # === Paso 7: Visualizaciones comparativas ===
    print(f"\n🔹 Paso {paso}: Generando gráficas comparativas...")
    comparativas_dir = results_dir / "comparativas"
    generar_visualizaciones(pre_csv=ora_pre_result, post_csv=ora_post_result, output_dir=comparativas_dir)
    print("📊 Gráficas comparativas generadas correctamente.")
    paso += 1

    # === Final ===
    print("\n🎯 Flujo de trabajo completado correctamente.")
    print(f"\n📁 Resultados principales:")
    print(f"   - ORA inicial:           {ora_pre_dir}")
    print(f"   - Propagación DIAMOnD:   {diamond_results_dir}")
    print(f"   - ORA posterior:         {ora_post_dir}")
    print(f"   - Comparativas:          {comparativas_dir}")


if __name__ == "__main__":
    main()