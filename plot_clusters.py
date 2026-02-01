import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Vizualizacija K-means klastera")

parser.add_argument(
    "--input",
    dest="input_C_file",
    type=str,
    default=None,
    help="Ulazni fajl za C program (opcion)"
)

parser.add_argument(
    "--output",
    dest="output_C_file",
    type=str,
    required=True,
    help="Izlazni fajl sa klaster podacima (obavezno)"
)

parser.add_argument(
    "--save",
    dest="save_file",
    type=str,
    required=True,
    help="Putanja gde se čuva slika (obavezno)"
)

args = parser.parse_args()

input_C_file = args.input_C_file
output_C_file = args.output_C_file
save_file = args.save_file

# Naziv fajla
# input_C_file = 'input_files/points_1000000.txt'

# output_C_file = 'output_files/clusters_data.txt'
# output_C_file = 'output_files/cluster_data_omp.txt'
# output_C_file = 'output_files/cluster_data_mpi.txt'
# output_C_file = 'output_files/clusters_500000_2.txt'

# save_file = 'images/cluster_data_seq_500000_2.png'
# save_file = 'images/cluster_data_omp.png'
# save_file = 'images/cluster_data_mpi.png'

# from ctypes import *
# so_file = "/home/klasick/faks/ubuntu-copy/semestar-1/HPC/projekat/k-means/sequential.so"
# my_functions = CDLL(so_file)

# if input_C_file is not None:
#     ret = my_functions.compute(
#         input_C_file.encode("utf-8"),
#         output_C_file.encode("utf-8")
#     )

clusters = []
current_cluster = []

with open(output_C_file, "r") as f:
    for line in f:
        line = line.strip()
        if not line:  # prazna linija znači kraj klastera
            if current_cluster:
                clusters.append(current_cluster)
                current_cluster = []
        elif line.startswith("Cluster"):
            continue  # preskoči liniju "Cluster X:"
        else:
            x_str, y_str = line.split()
            x, y = float(x_str), float(y_str)
            current_cluster.append((x, y))

# Dodaj poslednji klaster ako postoji
if current_cluster:
    clusters.append(current_cluster)

# Crtanje klastera
plt.figure(figsize=(8, 6))
colors = plt.cm.get_cmap('tab20', len(clusters))  # do 20 različitih boja

for idx, cluster in enumerate(clusters):
    xs, ys = zip(*cluster)  # raspakuj x i y koordinate
    plt.scatter(xs, ys, s=10, color=colors(idx), label=f"Cluster {idx}")

plt.xlabel("X")
plt.ylabel("Y")
plt.title("Vizualizacija klastera")
plt.legend()
plt.grid(True)
plt.savefig(save_file, dpi=1000)
print(f"Graf sačuvan u {save_file}")
# plt.show()
