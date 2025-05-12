import os
import time
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import sys

# Agregar carpeta del algoritmo original al path
sys.path.append(os.path.join(os.path.dirname(__file__), "parall_BFOA-main"))
from parallel_BFOA import run_bfoa  # Usa la función real

# Cargar archivo multifasta
fasta_file = "multifasta.fasta"
sequences = list(SeqIO.parse(fasta_file, "fasta"))
print(f"Se cargaron {len(sequences)} secuencias.")

# Ejecutar 30 corridas reales
resultados = []
for i in range(30):
    print(f"Corrida {i+1}/30...")
    start_time = time.time()
    fitness, iteraciones, blosum = run_bfoa()
    duracion = round(time.time() - start_time, 4)
    resultados.append({
        "Corrida": i + 1,
        "Fitness": fitness,
        "Tiempo (s)": duracion,
        "Iteraciones": iteraciones,
        "BLOSUM Score": blosum
    })

# Crear DataFrame y guardar CSV
df = pd.DataFrame(resultados)
df.to_csv("resultados_bfoa.csv", index=False)
print("✅ Se guardó el archivo resultados_bfoa.csv")

# Graficar resultados
plt.figure()
df["Fitness"].plot(kind="line", title="Fitness por Corrida")
plt.xlabel("Corrida")
plt.ylabel("Fitness")
plt.grid(True)
plt.savefig("grafica_fitness.png")

plt.figure()
df["Tiempo (s)"].plot(kind="bar", title="Tiempo de ejecución por Corrida")
plt.xlabel("Corrida")
plt.ylabel("Tiempo (s)")
plt.grid(True)
plt.tight_layout()
plt.savefig("grafica_tiempo.png")

plt.figure()
df["Iteraciones"].plot(kind="line", title="Iteraciones por Corrida", color="orange")
plt.xlabel("Corrida")
plt.ylabel("Iteraciones")
plt.grid(True)
plt.savefig("grafica_iteraciones.png")

plt.figure()
df["BLOSUM Score"].plot(kind="line", title="BLOSUM Score por Corrida", color="green")
plt.xlabel("Corrida")
plt.ylabel("BLOSUM Score")
plt.grid(True)
plt.savefig("grafica_blosum.png")

print("📈 Gráficas guardadas como PNG.")
