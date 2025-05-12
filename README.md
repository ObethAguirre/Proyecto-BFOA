# BFOA - Alineamiento de Secuencias con Algoritmo de Forrajeo Bacteriano

Este proyecto implementa y analiza un algoritmo de optimizacion por forrajeo bacteriano (BFOA) aplicado al alineamiento de secuencias usando la matriz BLOSUM62.

## Estructura

- `multifasta.fasta`: Archivo de entrada con las secuencias.
- `parall_BFOA-main/parallel_BFOA.py`: Algoritmo original.
- `analisis_bfoa.py`: Script que ejecuta 30 corridas del algoritmo y genera metricas y graficas.
- `resultados_bfoa.csv`: Datos generados con las 30 corridas.
- `reporte_graficas_bfoa.png`: Visualizacion de resultados.

## Requisitos

Instala las dependencias con pip:

```bash
pip install biopython matplotlib pandas python-docx
