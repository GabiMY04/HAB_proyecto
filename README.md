# 🌱 🧬 Análisis de propagación génica en *Arabidopsis thaliana*

Este proyecto implementa un flujo completo de análisis funcional y propagación en redes génicas en *Arabidopsis thaliana*, combinando el algoritmo **DIAMOnD (_Disease Module Detection_)** con un análisis funcional de sobrerrepresentación (**ORA**, _Over-Representation Analysis_).
El objetivo del flujo es identificar nuevos genes funcionalmente asociados a un conjunto inicial de genes semilla —provenientes de resultados de expresión diferencial— evaluando su conectividad dentro de la red y su relevancia biológica tras la propagación.

En este flujo se combinan dos niveles de análisis complementarios:
- El estructural, basado en la topología de la red mediante el algoritmo DIAMOnD, que permite expandir módulos génicos a partir de las semillas iniciales.
- El funcional, basado en el análisis de sobrerrepresentación (ORA), que evalúa los procesos biológicos enriquecidos antes y después de la propagación para revelar nuevas asociaciones funcionales entre los genes.

Este repositorio contiene un flujo completo en Python que integra una implementación
propia de la propagación DIAMOnD, el análisis funcional con `STRINGdb`, y la generación
de gráficas comparativas que permiten visualizar las diferencias funcionales pre y post-propagación.

---

## ⛓️ Descripción general del flujo

El proceso completo está completamente automatizado mediante el script principal `ejecutar_pipeline.py`, e incluye **siete pasos secuenciales**:

1. **Procesamiento de DEGs:** filtra genes diferencialmente expresados y genera listas separadas para genes sobre- y subexpresados.  
2. **Conversión a formato `STRINGdb`:** adapta los identificadores al prefijo taxonómico correspondiente a *Arabidopsis thaliana* (`3702.`).  
3. **Descarga (o reutilización) de red STRING:** obtiene o usa una red de interacciones proteína-proteína (PPI) centrada en los genes semilla.  
4. **Análisis funcional inicial (ORA):** evalúa las funciones biológicas enriquecidas entre los genes semilla.  
5. **Propagación con DIAMOnD:** expande el conjunto inicial añadiendo genes significativamente conectados.  
6. **Análisis funcional posterior (ORA):** repite el análisis funcional sobre los genes propagados.  
7. **Comparación visual pre/post:** genera gráficas y diagramas que muestran los cambios en la significancia funcional tras la propagación.

> 💡 El flujo es completamente modular: cada paso puede ejecutarse y evaluarse de forma independiente desde el script correspondiente.

---

## 📁 Estructura del repositorio

```
/analisis-arabidopsis/
├── data/
│   └── Allcontrasts_GLM-Treat_*.tsv      # Archivo de entrada (expresión diferencial)
│
├── scripts/
│   ├── procesar_DEGs.py                  # Filtra genes diferencialmente expresados
│   ├── convertir_ids_string.py           # Convierte AGI → STRINGdb
│   ├── descargar_red_string.py           # Descarga la red de interacciones
│   ├── analisis_funcional.py             # Implementación propia del ORA mediante STRINGdb
│   ├── diamond.py                        # Implementación propia del algoritmo DIAMOnD  
│   └── comparar_enriquecimientos.py      # Visualización comparativa de resultados  
│
├── results/
│   ├── ORA_semillas/                     # Resultados del ORA inicial
│   ├── diamond_propagation/              # Resultados de la propagación DIAMOnD
│   ├── ORA_diamond/                      # Resultados del ORA posterior
│   └── comparativas/                     # Gráficas comparativas
│    
│
├── ejecutar_pipeline.py                  # Script principal (flujo completo)
├── README.md
└── requirements.txt
```

---

## 🔍 Implementación de ORA

El ORA clásico emplea la prueba exacta de Fisher, una herramienta estadística que evalúa
si existe una asociación significativa entre dos variables categóricas —en este caso, los genes de interés y las categorías
funcionales a las que pertenecen—. En este flujo, la implementación del ORA utiliza precisamente este enfoque estadístico,
ejecutado a través de la API oficial de `STRINGdb`, que aplica internamente la prueba de Fisher para evaluar el enriquecimiento
de cada categoría funcional.  

El análisis se restringe exclusivamente a los términos de la base de datos Gene Ontology (GO), lo que permite obtener una descripción
funcional clara y estandarizada de los genes analizados. Los resultados se corrigen mediante el FDR (_False Discovery Rate_),
una versión ajustada del *p*-valor que controla la proporción esperada de falsos positivos, garantizando así
una interpretación estadística más robusta y confiable.

Finalmente, los identificadores GO (por ejemplo, `GO:0015979`) se traducen automáticamente a sus nombres descriptivos
mediante la ontología oficial de Gene Ontology (`go-basic.obo`), produciendo una tabla interpretable de las funciones
biológicas correspondientes.. El análisis genera dos salidas:

- `enrichment_results.csv`: tabla con las categorías GO, sus valores FDR y genes asociados.  
- `enrichment_plot.png`: gráfico de barras con las categorías más significativamente enriquecidas.

## ♦️ Implementación de DIAMOnD

El método DIAMOnD considera que los genes relacionados con una misma enfermedad o proceso biológico tienden a agruparse en módulos densamente
conectado. La implementación del algoritmo fue desarrollada de forma personalizada,
adaptada a la estructura y objetivos del flujo de trabajo. El método expande iterativamente un conjunto de genes semilla dentro de una red génica,
añadiendo en cada iteración el gen más significativamente conectado según una prueba hipergeométrica.
La implementación es completamente determinista: para un mismo conjunto de semillas y red de entrada, el resultado del módulo generado será siempre idéntico.
Esto garantiza la reproducibilidad completa de los resultados.

Durante la ejecución, el algoritmo se detiene automáticamente cuando se cumple
alguna de las siguientes condiciones:
 - 1.	No existen más candidatos posibles.
 - 2.	El mejor _p_-valor deja de ser estadísticamente significativo (_p_ > 0.05).
 - 3.	Se alcanzan 100 genes añadidos al módulo.
- 4.	Se completan 200 iteraciones del proceso.

Estos criterios garantizan una expansión controlada del módulo y evitan la incorporación de nodos con baja relevancia estadística.
Los parámetros internos del modelo se encuentran fijados y no deben modificarse,
con el fin de mantener la reproducibilidad de los resultados.

Cada ejecución de DIAMOnD genera dos archivos de salida en el directorio `results/diamond_propagation/`:
- `diamond_results.csv`: tabla con los genes del módulo (semillas + añadidos), sus _p_-valores y número de conexiones.  
- `diamond_genes.txt`: lista de genes del módulo completo (uno por línea).

---

## 🚀 Manual de uso

Clonar el repositorio y ejecutar el pipeline completo:

```bash
git clone https://github.com/GabiMY04/HAB_proyecto analisis-arabidopsis
cd analisis-arabidopsis
pip install -r requirements.txt
python ejecutar_pipeline.py
```

Parámetros de ejecución:

| Parámetro     | Descripción                                                                                                 | Opcional | Valor por defecto |
|---------------|-------------------------------------------------------------------------------------------------------------|-----------|-------------------|
| `--input_dir` | Directorio que contiene los archivos de entrada (`Allcontrasts_GLM-Treat_*.tsv`). | ✅        | `data/`           |
| `--output_dir` | Directorio raíz donde se guardarán todos los resultados generados por el flujo.                             | ✅        | `results/`        |

> 💡 Los archivos de entrada deben ubicarse en la carpeta data/.
Si se desea usar otro resultado de expresión diferencial, basta con reemplazarlo antes de ejecutar el script.


## 📊 Resultados

Los resultados se guardan automáticamente en el directorio especificado mediante el parámetro `--output_dir` (por defecto, `results/`) y se organizan por tipo de análisis:

- `ORA_semillas/` y `ORA_diamond/`: contienen los resultados del análisis funcional (ORA) realizado antes y después de la propagación, respectivamente.  
  Cada carpeta incluye:  
  - `enrichment_results.csv`: tabla con las categorías funcionales enriquecidas, sus valores FDR y los genes asociados.  
  - `enrichment_plot.png`: gráfico de barras con las categorías funcionales más significativas obtenidas en el análisis.  

- `diamond_propagation/`: almacena los resultados del algoritmo DIAMOnD, incluyendo los genes añadidos al módulo y sus respectivos _p_-valores (`diamond_results.csv`), así como la lista completa de genes propagados (`diamond_genes.txt`).  

- `comparativas/`: contiene las visualizaciones que comparan los resultados del análisis funcional antes y después de aplicar DIAMOnD:
  - `ora_barplot_delta.png`: gráfico de barras comparativo que muestra, para las categorías compartidas entre ambos análisis, la diferencia en
  significancia estadística (_-log10 FDR) entre el ORA pre y post-propagación. Permite observar qué funciones ganan o pierden relevancia tras el DIAMOnD.  
  - `ora_barplot_nuevas.png`: gráfico de barras con las nuevas categorías funcionales que aparecen únicamente
  después de la propagación, junto con su nivel de significancia.  
  - `ora_venn.png`: diagrama de Venn que representa el grado de solapamiento entre los términos funcionales enriquecidos
  antes y después de la propagación, distinguiendo las categorías compartidas y las exclusivas de cada análisis.  

## ⚙️ Dependencias

Las librerías necesarias para ejecutar el análisis se encuentran en `requirements.txt`:

```
networkx
pandas
numpy
scipy
matplotlib
goatools
stringdb
matplotlib-venn
```

## 📚 Referencias

- Ghiassian et al., *A DIseAse MOdule Detection (DIAMOnD) Algorithm derived from a Systematic Analysis of Connectivity Patterns in Disease Networks.* **PLOS Computational Biology**, 2015.  
- Szklarczyk et al., *STRING v12: Protein–protein association networks in 2024.* **Nucleic Acids Research**, 2024.  
- The Gene Ontology Consortium, *Gene Ontology: the unified resource for gene annotation.* **Nature Genetics**, 2021.  