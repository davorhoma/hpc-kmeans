## Pregled projekta

K-means je algoritam za klasterovanje koji deli skup podataka u K klastera. Svaka tačka se dodeljuje najbližem centroidu, a centroidi se ažuriraju iterativno dok ne konvergiraju.

Ovaj projekat sadrži tri implementacije:

1.  **Sequential** (`sequential.c`) - bazna sekvencijalna verzija
2.  **OpenMP** (`omp_parallel.c`) - paralelizovana verzija korišćenjem OpenMP (paralelizacija na nivou niti)
3.  **MPI** (`mpi_parallel.c`) - distribuirana verzija korišćenjem MPI (paralelizacija na nivou procesa)

## Preduslov

Potrebni alati i biblioteke:

```bash
# GCC kompajler
sudo apt-get install build-essential

# OpenMP (obično dolazi sa GCC-om)
# Proveriti: gcc -fopenmp --version

# OpenMPI
sudo apt-get install libopenmpi-dev openmpi-common openmpi-bin

# bc za testiranje koristeći skriptu
sudo apt-get install bc

```

## Kompilacija
```bash
# Sekvencijalna verzija
gcc -o seq sequential.c -lm -fopenmp

# OpenMP verzija
gcc -o omp omp_parallel.c -fopenmp -lm

# MPI verzija
mpicc -o mpi mpi_parallel.c -fopenmp -lm

```

## Pokretanje

Svi programi primaju iste argumente:

```bash
<program> <input_file> <output_file> <K>

```

-   `input_file` - putanja do fajla sa tačkama
-   `output_file` - putanja do izlaznog fajla gde će biti sačuvani klasteri
-   `K` - broj klastera

### Sekvencijalna verzija

```bash
./seq input_files/points.txt output_files/clusters.txt 4

```

### OpenMP verzija

```bash
# Podrazumevano koristi sve dostupne jezgre
./omp input_files/points.txt output_files/clusters_omp.txt 4

# Specificiranje broja niti
export OMP_NUM_THREADS=4
./omp input_files/points.txt output_files/clusters_omp.txt 4

```

### MPI verzija

```bash
# Pokretanje sa 4 procesa
mpirun -np 4 ./mpi input_files/points.txt output_files/clusters_mpi.txt 4

# Na sistemima sa hiperthreading-om
mpirun --use-hwthread-cpus -np 8 ./mpi input_files/points.txt output_files/clusters_mpi.txt 4

```

## Testiranje performansi

Skripta `performance_test.sh` automatizuje testiranje performansi za različite veličine ulaznih podataka i različit broj niti/procesa.

### Priprema testnih podataka

Prvo kreirajte direktorijum i generišite test fajlove:

```bash
mkdir -p input_files output_files

# Generisanje test fajlova različitih veličina (opciono - možete koristiti svoje podatke)
# points_400.txt, points_4000.txt, points_40000.txt, points_400000.txt, points_4000000.txt

```

### Konfiguracija test skripte

Podesite parametre u `performance_test.sh`:

```bash
INPUT_FILES=(
    "input_files/points_400.txt"
    "input_files/points_4000.txt"
    "input_files/points_40000.txt"
    "input_files/points_400000.txt"
    "input_files/points_4000000.txt"
)

THREADS=(2 4 8 16)  # Broj niti/procesa za testiranje
K=4                  # Broj klastera
REPEAT=3             # Broj ponavljanja za svaki test

```

### Pokretanje testova

```bash
# Dajte skripti izvršna prava
chmod +x performance_test.sh

# Pokrenite testiranje
./performance_test.sh

```

Rezultati će biti sačuvani u `results_scaling.csv` fajlu u sledećem formatu:

```csv
Implementation,InputSize,Workers,BestTime
Sequential,400,1,0.002345
OpenMP,400,2,0.001567
OpenMP,400,4,0.001234
MPI,400,2,0.001890
...

```

## Struktura projekta

```
.
├── sequential.c                # Sekvencijalna implementacija
├── omp_parallel.c              # OpenMP implementacija
├── mpi_parallel.c              # OpenMPI implementacija
├── performance_test.sh         # Skripta za testiranje performansi
├── generator.py                # Skripta za generisanje test podataka
├── plot_clusters.py            # Skripta za vizualizaciju klastera
├── plot_scaling.py             # Skripta za analizu performansi
├── input_files/                # Direktorijum za ulazne fajlove
│   ├── points_400.txt
│   └── ...
├── output_files/               # Direktorijum za izlazne fajlove
│   ├── clusters.txt
│   └── ...
└── images/                     # Direktorijum za grafike
    ├── scaling_omp.png
    ├── scaling_mpi.png
    └── ...

```

## Format ulaznih podataka

Ulazni fajlovi mogu biti u `.txt` ili `.csv` formatu sa koordinatama tačaka:

```
x1,y1
x2,y2
x3,y3
...

```

Primer (`points.txt`):

```
0.234,0.567
0.891,0.123
0.456,0.789
0.112,0.445

```

## Izlazni fajlovi

Svaki program generiše izlazni fajl sa formatom:

```
Cluster 0:
0.234000 0.567000
0.456000 0.789000

Cluster 1:
0.891000 0.123000
0.112000 0.445000

...

```

## Python skripte

Projekat uključuje tri Python skripte za generisanje podataka i vizualizaciju rezultata.

### 1. Generisanje testnih podataka (`generator.py`)

Skripta za generisanje sintetičkih podataka sa kontrolisanim brojem klastera.

**Upotreba:**

```bash
# Osnovni primer - generiše 3 klastera sa po 100 tačaka
python generator.py --out input_files/points.csv

# Napredni primer
python generator.py \
    --clusters 5 \
    --points 10000 \
    --spread 0.3 \
    --distance 8 \
    --out input_files/points_50000.txt \
    --txt \
    --plot

```

**Parametri:**

-   `--clusters` - broj klastera (default: 3)
-   `--points` - broj tačaka po klasteru (default: 100)
-   `--spread` - raspršenost tačaka unutar klastera (default: 0.5)
-   `--distance` - rastojanje između centroida klastera (default: 5)
-   `--out` - izlazni fajl (default: points.csv)
-   `--txt` - sačuvaj kao .txt umesto .csv
-   `--plot` - prikaži i sačuvaj grafik (generator_clusters.png)

**Primer generisanja svih test fajlova:**

```bash
python generator.py --clusters 4 --points 100 --out input_files/points_400.txt --txt
python generator.py --clusters 4 --points 1000 --out input_files/points_4000.txt --txt
python generator.py --clusters 4 --points 10000 --out input_files/points_40000.txt --txt
python generator.py --clusters 4 --points 100000 --out input_files/points_400000.txt --txt
python generator.py --clusters 4 --points 1000000 --out input_files/points_4000000.txt --txt

```

### 2. Vizualizacija klastera (`plot_clusters.py`)

Skripta za vizualizaciju rezultata klasterovanja.

**Upotreba:**

```bash
python plot_clusters.py \
    --output output_files/clusters.txt \
    --save images/clusters_visualization.png
```

**Parametri (obavezni):**
-   `--output` - fajl sa rezultatima klasterovanja
-   `--save` - putanja gde se čuva slika

Skripta čita izlazni fajl i kreira scatter plot gde svaki klaster ima različitu boju.

### 3. Analiza performansi (`plot_scaling.py`)

Skripta za vizualizaciju scaling performansi na osnovu `results_scaling.csv`.

**Upotreba:**

```bash
# Napravite direktorijum za slike
mkdir -p images

# Pokrenite analizu
python plot_scaling.py

```

Skripta automatski generiše **tri grafikona**:

1.  **`images/scaling_omp.png`** - OpenMP scaling u poređenju sa sekvencijalnom verzijom
2.  **`images/scaling_mpi.png`** - MPI scaling u poređenju sa sekvencijalnom verzijom
3.  **`images/scaling_omp_vs_mpi.png`** - Direktno poređenje OpenMP i MPI performansi

Grafici koriste logaritamsku skalu za obe ose kako bi pokazali scaling karakteristike za različite veličine dataset-a.

**Napomena:** Skripta očekuje da postoji `results_scaling.csv` fajl generisan sa `performance_test.sh`.

## Kompletna sekvenca testiranja

Evo kompletnog primera kako da generišete podatke, pokrenete testove i vizualizujete rezultate:

```bash
# 1. Kreirajte direktorijume
mkdir -p input_files output_files images

# 2. Generišite test podatke
python generator.py --clusters 4 --points 100 --out input_files/points_400.txt --txt
python generator.py --clusters 4 --points 1000 --out input_files/points_4000.txt --txt
python generator.py --clusters 4 --points 10000 --out input_files/points_40000.txt --txt
python generator.py --clusters 4 --points 100000 --out input_files/points_400000.txt --txt
python generator.py --clusters 4 --points 1000000 --out input_files/points_4000000.txt --txt

# 3. Kompajlirajte programe
gcc -o seq sequential.c -lm -fopenmp
gcc -o omp omp_parallel.c -fopenmp -lm
mpicc -o mpi mpi_parallel.c -fopenmp -lm

# 4. Pokrenite performance testove
chmod +x performance_test.sh
./performance_test.sh

# 5. Vizualizujte rezultate
python plot_scaling.py

# 6. Vizualizujte klastering na jednom primeru
./seq input_files/points_4000.txt output_files/clusters_seq.txt 4
python plot_clusters.py --output output_files/clusters_seq.txt --save images/clusters_example.png

```