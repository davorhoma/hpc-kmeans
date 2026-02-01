import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ---- Ucitaj CSV ----
df = pd.read_csv("results_scaling.csv")

# ---- Dataset-ovi koje hocemo da prikazemo ----
datasets = [400, 4000, 40000, 400000, 4000000]
omp_save_file = "images/scaling_omp.png"
mpi_save_file = "images/scaling_mpi.png"

# ---- X osa: log2 scale ----
workers = [1, 2, 4, 8, 16]
x = np.log2(workers)

# ---- 1️⃣ Sekvencijalni + OpenMP ----
plt.figure(figsize=(10, 6))
for n in datasets:
    subset = df[(df['InputSize']==n) & (df['Implementation'].isin(['Sequential','OpenMP']))]
    # Nadji sekvencijalno vreme (Workers=1)
    seq_time = subset[subset['Implementation']=='Sequential']['BestTime'].values[0]
    # OpenMP vreme po broju niti
    omp_times = []
    for w in workers:
        row = subset[(subset['Implementation']=='OpenMP') & (subset['Workers']==w)]
        if len(row)==0:
            # Ako nema podatka (npr. mali dataset sa vise niti)
            omp_times.append(np.nan)
        else:
            omp_times.append(row['BestTime'].values[0])
    # Kombinujemo X i Y
    times = [seq_time] + omp_times[1:]  # Za 1 thread uzimamo sekvencijalno vreme
    plt.plot(x, times, marker='o', label=f"{n} points")

plt.xticks(x, workers)
plt.xlabel("Number of threads")
plt.ylabel("Execution time (s)")
plt.yscale("log", base=10)
plt.title("OpenMP Scaling vs Sequential")
plt.grid(True, which="both", ls="--")
plt.legend()
plt.tight_layout()
plt.savefig(omp_save_file, dpi=1000)
print(f"Graf sačuvan u {omp_save_file}")

# ---- 2️⃣ Sekvencijalni + MPI ----
plt.figure(figsize=(10, 6))
for n in datasets:
    subset = df[(df['InputSize']==n) & (df['Implementation'].isin(['Sequential','MPI']))]
    # Nadji sekvencijalno vreme (Workers=1)
    seq_time = subset[subset['Implementation']=='Sequential']['BestTime'].values[0]
    # MPI vreme po broju procesa
    mpi_times = []
    for w in workers:
        row = subset[(subset['Implementation']=='MPI') & (subset['Workers']==w)]
        if len(row)==0:
            mpi_times.append(np.nan)
        else:
            mpi_times.append(row['BestTime'].values[0])
    times = [seq_time] + mpi_times[1:]
    plt.plot(x, times, marker='o', label=f"{n} points")

plt.xticks(x, workers)
plt.xlabel("Number of processes")
plt.ylabel("Execution time (s)")
plt.yscale("log", base=10)
plt.title("MPI Scaling vs Sequential")
plt.grid(True, which="both", ls="--")
plt.legend()
plt.tight_layout()
plt.savefig(mpi_save_file, dpi=1000)
print(f"Graf sačuvan u {mpi_save_file}")

# ---- 3️⃣ Poređenje OpenMP vs MPI ----

compare_save_file = "images/scaling_omp_vs_mpi.png"

plt.figure(figsize=(11, 7))

colors = ['blue', 'orange', 'green', 'red', 'purple']
for n in datasets:
    subset = df[df['InputSize'] == n]

    # Sekvencijalno vreme
    seq_time = subset[subset['Implementation']=='Sequential']['BestTime'].values[0]

    omp_times = []
    mpi_times = []

    for w in workers:
        # OpenMP
        omp_row = subset[(subset['Implementation']=='OpenMP') & (subset['Workers']==w)]
        if len(omp_row)==0:
            omp_times.append(np.nan)
        else:
            omp_times.append(omp_row['BestTime'].values[0])

        # MPI
        mpi_row = subset[(subset['Implementation']=='MPI') & (subset['Workers']==w)]
        if len(mpi_row)==0:
            mpi_times.append(np.nan)
        else:
            mpi_times.append(mpi_row['BestTime'].values[0])

    # Ubacujemo sekvencijalno vreme za 1 worker
    omp_plot = [seq_time] + omp_times[1:]
    mpi_plot = [seq_time] + mpi_times[1:]

    # Crtanje
    plt.plot(x, omp_plot, marker='o', linestyle='-', color=colors[datasets.index(n)], label=f"OpenMP {n}")
    plt.plot(x, mpi_plot, marker='s', linestyle='--', color=colors[datasets.index(n)], label=f"MPI {n}")

plt.xticks(x, workers)
plt.xlabel("Number of threads / processes")
plt.ylabel("Execution time (s)")
plt.yscale("log", base=10)
plt.title("OpenMP vs MPI Scaling Comparison")
plt.grid(True, which="both", ls="--")
plt.legend(ncol=2)
plt.tight_layout()
plt.savefig(compare_save_file, dpi=1000)

print(f"Poređenje OpenMP vs MPI sačuvano u {compare_save_file}")
