#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

struct Point {
    double x;
    double y;
};

struct Cluster {
    struct Point* points;
    int size;
    int capacity;
};

struct ClusterAccum {
    double sum_x;
    double sum_y;
    int count;
};

MPI_Datatype MPI_POINT;
MPI_Datatype MPI_CLUSTER_ACCUM;
MPI_Op MPI_CLUSTER_SUM;

void create_mpi_cluster_accum_type();
void create_mpi_point_type();
void cluster_accum_sum(void*, void*, int*, MPI_Datatype*);

double euclidean_distance(struct Point*, struct Point*);
int* kmeans_mpi(int, struct Point*, int, struct Point*);

int load_points_from_file(const char*, struct Point**);
void save_clusters_to_file(const char*, struct Cluster*, int K);

int compute(const char* input_filename, const char* output_filename, int K, int rank, int size) {
    srand(time(NULL) + rank);

    struct Point* all_points = NULL;
    struct Point* local_points = NULL;

    int n = 0;
    int local_n = 0;

    int* sendcounts = NULL;
    int* displs = NULL;

    if (rank == 0) {
        n = load_points_from_file(input_filename, &all_points);
        if (n <= 0) {
            printf("No points loaded\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        sendcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));

        int base = n / size;
        int rem = n % size;

        for (int i = 0; i < size; i++) {
            sendcounts[i] = base + (i < rem ? 1 : 0);
            displs[i] = (i == 0) ? 0 : displs[i - 1] + sendcounts[i - 1];
        }
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    local_n = n / size;
    if (rank < (n % size))
        local_n++;

    local_points = malloc(local_n * sizeof(struct Point));
    if (!local_points) {
        perror("malloc failed");
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    MPI_Scatterv(all_points, sendcounts, displs, MPI_POINT, local_points, local_n, MPI_POINT, 0, MPI_COMM_WORLD);

    double start, end;
    if (rank == 0) start = MPI_Wtime();

    struct Point centroids[K];
    if (rank == 0) {
        for (int i = 0; i < K; i++) {
            centroids[i].x = (double)rand() / RAND_MAX;
            centroids[i].y = (double)rand() / RAND_MAX;
        }
    }

    MPI_Bcast(centroids, K, MPI_POINT, 0, MPI_COMM_WORLD);

    int* labels = kmeans_mpi(K, local_points, local_n, centroids);

    if (rank == 0) {
        end = MPI_Wtime();
        printf("Time elapsed: %lf\n", end - start);
    }

    int local_counts[K];
    memset(local_counts, 0, K * sizeof(int));
    for (int i = 0; i < local_n; ++i) {
        local_counts[labels[i]]++;
    }

    int global_counts[K];
    MPI_Reduce(local_counts, global_counts, K, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    struct Cluster* clusters = NULL;

    if (rank == 0) {
        clusters = malloc(K * sizeof(struct Cluster));

        for (int k = 0; k < K; ++k) {
            clusters[k].size = global_counts[k];
            clusters[k].capacity = global_counts[k];
            clusters[k].points = malloc(global_counts[k] * sizeof(struct Point));
        }
    }

    for (int k = 0; k < K; ++k) {
        struct Point* local_buf = malloc(local_counts[k] * sizeof(struct Point));
        int idx = 0;

        for (int i = 0; i < local_n; ++i) {
            if (labels[i] == k) {
                local_buf[idx++] = local_points[i];
            }
        }

        int* recv_counts = NULL;
        int* rdispls = NULL;

        if (rank == 0) {
            recv_counts = malloc(size * sizeof(int));
        }

        MPI_Gather(&local_counts[k], 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            rdispls = malloc(size * sizeof(int));
            rdispls[0] = 0;
            for (int i = 1; i < size; ++i) {
                rdispls[i] = rdispls[i - 1] + recv_counts[i - 1];
            }
        }

        MPI_Gatherv(local_buf, local_counts[k], MPI_POINT, rank == 0 ? clusters[k].points : NULL, recv_counts, rdispls, MPI_POINT, 0, MPI_COMM_WORLD);

        free(local_buf);
        if (rank == 0) {
            free(recv_counts);
            free(rdispls);
        }
    }

    if (rank == 0) {
        save_clusters_to_file(output_filename, clusters, K);

        for (int i = 0; i < K; ++i) {
            printf("Centroid %d: (%.4f, %.4f)\n",
                   i, centroids[i].x, centroids[i].y);
        }
    }

    free(labels);
    free(local_points);

    if (rank == 0) {
        for (int k = 0; k < K; ++k)
            free(clusters[k].points);

        free(clusters);
        free(all_points);
        free(sendcounts);
        free(displs);
    }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        printf("Usage: %s <input_file> <output_file> <K>\n", argv[0]);
        printf("Example: %s input_files/points.txt output_files/clusters.txt 4\n", argv[0]);
        return 1;
    }

    const char* input_file = argv[1];
    const char* output_file = argv[2];
    int K = atoi(argv[3]);

    if (K <= 0) {
        printf("Error: K must be positive\n");
        return 1;
    }

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        printf("Running with %d MPI processes\n", size);
    }

    create_mpi_point_type();
    create_mpi_cluster_accum_type();
    MPI_Op_create(cluster_accum_sum, 1, &MPI_CLUSTER_SUM);

    compute(input_file, output_file, K, rank, size);

    MPI_Op_free(&MPI_CLUSTER_SUM);
    MPI_Type_free(&MPI_CLUSTER_ACCUM);
    MPI_Type_free(&MPI_POINT);
    MPI_Finalize();
    return 0;
}

double euclidean_distance(struct Point* a, struct Point* b) {
    double dx = a->x - b->x;
    double dy = a->y - b->y;
    return sqrt(dx * dx + dy * dy);
}

void create_mpi_cluster_accum_type() {
    int block_lengths[3] = {1, 1, 1};
    MPI_Aint offsets[3];
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};

    offsets[0] = offsetof(struct ClusterAccum, sum_x);
    offsets[1] = offsetof(struct ClusterAccum, sum_y);
    offsets[2] = offsetof(struct ClusterAccum, count);

    MPI_Type_create_struct(3, block_lengths, offsets, types, &MPI_CLUSTER_ACCUM);
    MPI_Type_commit(&MPI_CLUSTER_ACCUM);
}

void create_mpi_point_type() {
    int block_lengths[2] = {1, 1};
    MPI_Aint offsets[2];
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};

    offsets[0] = offsetof(struct Point, x);
    offsets[1] = offsetof(struct Point, y);

    MPI_Type_create_struct(2, block_lengths, offsets, types, &MPI_POINT);
    MPI_Type_commit(&MPI_POINT);
}

void cluster_accum_sum(void* in, void* inout, int* len, MPI_Datatype* dptr) {
    struct ClusterAccum* a = in;
    struct ClusterAccum* b = inout;

    for (int i = 0; i < *len; i++) {
        b[i].sum_x += a[i].sum_x;
        b[i].sum_y += a[i].sum_y;
        b[i].count += a[i].count;
    }
}

int* kmeans_mpi(int K, struct Point* local_points, int local_n, struct Point* centroids) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    struct ClusterAccum local[K];
    struct ClusterAccum global[K];
    int* labels = malloc(local_n * sizeof(int));

    int max_iters = 100;
    double tolerance = 1e-4;

    for (int iter = 0; iter < max_iters; ++iter) {
        for (int k = 0; k < K; ++k) {
            local[k].sum_x = 0.0;
            local[k].sum_y = 0.0;
            local[k].count = 0;
        }

        for (int i = 0; i < local_n; ++i) {
            int best = 0;
            double min_dist = euclidean_distance(&local_points[i], &centroids[0]);

            for (int k = 1; k < K; ++k) {
                double distance = euclidean_distance(&local_points[i], &centroids[k]);
                if (distance < min_dist) {
                    min_dist = distance;
                    best = k;
                }
            }

            local[best].sum_x += local_points[i].x;
            local[best].sum_y += local_points[i].y;
            local[best].count++;

            labels[i] = best;
        }

        MPI_Allreduce(local, global, K, MPI_CLUSTER_ACCUM, MPI_CLUSTER_SUM, MPI_COMM_WORLD);

        double max_shift = 0.0;

        for (int k = 0; k < K; ++k) {
            if (global[k].count > 0) {
                struct Point new_centroid;
                new_centroid.x = global[k].sum_x / global[k].count;
                new_centroid.y = global[k].sum_y / global[k].count;

                double shift = euclidean_distance(&centroids[k], &new_centroid);

                if (shift > max_shift) max_shift = shift;

                centroids[k] = new_centroid;
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, &max_shift, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if (max_shift < tolerance) break;
    }

    return labels;
}

int load_points_from_file(const char* filename, struct Point** points) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        printf("%s\n", filename);
        perror("Cannot open file");
        return -1;
    }

    int capacity = 100;
    int size = 0;
    *points = malloc(capacity * sizeof(struct Point));
    if (!*points) exit(EXIT_FAILURE);

    double x, y;
    while (fscanf(fp, "%lf,%lf", &x, &y) == 2) {
        if (size == capacity) {
            capacity *= 2;
            *points = realloc(*points, capacity * sizeof(struct Point));
            if (!*points) exit(EXIT_FAILURE);
        }

        (*points)[size].x = x;
        (*points)[size].y = y;
        size++;
    }

    fclose(fp);
    return size;
}

void save_clusters_to_file(const char* filename, struct Cluster* clusters, int K) {
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        perror("Cannot open output file");
        return;
    }

    for (int i = 0; i < K; ++i) {
        fprintf(fp, "Cluster %d:\n", i);
        for (int j = 0; j < clusters[i].size; j++) {
            fprintf(fp, "%.6f %.6f\n",
                    clusters[i].points[j].x,
                    clusters[i].points[j].y);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}
