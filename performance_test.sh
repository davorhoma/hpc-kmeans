#!/bin/bash

# ---------------- CONFIG ----------------

INPUT_FILES=(
    "input_files/points_400.txt"
    "input_files/points_4000.txt"
    "input_files/points_40000.txt"
    "input_files/points_400000.txt"
    "input_files/points_4000000.txt"
)

THREADS=(2 4 8 16)

OUTPUT_DIR="output_files"

K=4

SEQ_EXEC="./seq"
OMP_EXEC="./omp"
MPI_EXEC="./mpi"

REPEAT=3

RESULT_FILE="results_scaling.csv"

# ----------------------------------------

mkdir -p $OUTPUT_DIR

echo "Implementation,InputSize,Workers,BestTime" > $RESULT_FILE

extract_time () {
    echo "$1" | grep "Time elapsed" | awk '{print $3}'
}

get_best_time () {

    CMD="$1"

    BEST=9999999

    for ((i=1;i<=REPEAT;i++)); do

        OUTPUT=$(eval $CMD)

        TIME=$(extract_time "$OUTPUT")

        # echo "    Run $i: $TIME s"

        COMP=$(echo "$TIME < $BEST" | bc)

        if [ "$COMP" -eq 1 ]; then
            BEST=$TIME
        fi
    done

    echo $BEST
}

echo "======= K-MEANS SCALING BENCHMARK ======="

# ---------------- MAIN LOOP ----------------

for INPUT in "${INPUT_FILES[@]}"; do

    SIZE=$(echo $INPUT | grep -oE '[0-9]+')

    echo ""
    echo "======================================="
    echo "Dataset size: $SIZE"
    echo "======================================="

    # ---------- SEQUENTIAL ----------

    echo "Sequential:"

    CMD_SEQ="$SEQ_EXEC $INPUT $OUTPUT_DIR/seq_${SIZE}.txt $K"

    BEST_SEQ=$(get_best_time "$CMD_SEQ")

    echo "Best Sequential time: $BEST_SEQ s"

    echo "Sequential,$SIZE,1,$BEST_SEQ" >> $RESULT_FILE


    # ---------- OPENMP ----------

    echo ""
    echo "OpenMP scaling:"

    for T in "${THREADS[@]}"; do

        echo "  Threads = $T"

        export OMP_NUM_THREADS=$T

        CMD_OMP="$OMP_EXEC $INPUT $OUTPUT_DIR/omp_${SIZE}_${T}.txt $K"

        BEST_OMP=$(get_best_time "$CMD_OMP")

        echo "  Best OpenMP ($T): $BEST_OMP s"

        echo "OpenMP,$SIZE,$T,$BEST_OMP" >> $RESULT_FILE

    done


    # ---------- MPI ----------

    echo ""
    echo "MPI scaling:"

    for P in "${THREADS[@]}"; do

        echo "  Processes = $P"

        CMD_MPI="mpirun --use-hwthread-cpus -np $P $MPI_EXEC $INPUT $OUTPUT_DIR/mpi_${SIZE}_${P}.txt $K"

        BEST_MPI=$(get_best_time "$CMD_MPI")

        echo "  Best MPI ($P): $BEST_MPI s"

        echo "MPI,$SIZE,$P,$BEST_MPI" >> $RESULT_FILE

    done

done

echo ""
echo "======================================="
echo "All results saved to: $RESULT_FILE"
