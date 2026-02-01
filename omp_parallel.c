#include <math.h>
#include <omp.h>
#include <stdbool.h>
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

void add_point_to_cluster(struct Cluster*, struct Point);
bool are_points_equal(struct Point*, struct Point*);
double euclidean_distance(struct Point*, struct Point*);
struct Cluster* kmeans(int, struct Point[], int);
struct Point calculate_centroid(struct Point*, int, struct Point);

int load_points_from_file(const char*, struct Point**);
void save_clusters_to_file(const char*, struct Cluster*, int K);

int compute(const char* input_filename, const char* output_filename, int K) {
    srand(time(NULL));

    struct Point* points = NULL;
    int n = load_points_from_file(input_filename, &points);
    if (n <= 0) {
        printf("No points loaded\n");
        return 1;
    }

    // clock_t start = clock();
    double start = omp_get_wtime();
    struct Cluster* clusters = kmeans(K, points, n);
    // clock_t end = clock();
    double end = omp_get_wtime();

    // double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time elapsed: %lf\n", end - start);

    save_clusters_to_file(output_filename, clusters, K);

    for (int i = 0; i < K; i++) {
        free(clusters[i].points);
    }
    free(clusters);
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
    
    compute(input_file, output_file, K);
}

void add_point_to_cluster(struct Cluster* cluster, struct Point p) {
    if (cluster->size >= cluster->capacity) {
        cluster->capacity = (cluster->capacity == 0) ? 1 : cluster->capacity * 2;
        cluster->points = realloc(cluster->points, cluster->capacity * sizeof(struct Point));

        if (!cluster->points) {
            perror("realoc failed");
            exit(EXIT_FAILURE);
        }
    }

    cluster->points[cluster->size++] = p;
}

double euclidean_distance(struct Point* a, struct Point* b) {
    double dx = a->x - b->x;
    double dy = a->y - b->y;
    return sqrt(dx * dx + dy * dy);
}

struct Point calculate_centroid(struct Point* points, int size, struct Point old_centroid) {
    if (size == 0) return old_centroid;

    struct Point mean = {0.0, 0.0};

    for (int i = 0; i < size; i++) {
        mean.x += points[i].x;
        mean.y += points[i].y;
    }

    mean.x /= size;
    mean.y /= size;

    return mean;
}

struct Cluster* kmeans(int K, struct Point points[], int n) {
    struct Point centroids[K];
    struct Cluster* clusters = malloc(K * sizeof(struct Cluster));
    struct Point newCentroids[K];

    for (int i = 0; i < K; ++i) {
        clusters[i].points = malloc(n * sizeof(struct Point));
        clusters[i].size = 0;
        clusters[i].capacity = 0;

        centroids[i].x = (double)rand() / RAND_MAX;
        centroids[i].y = (double)rand() / RAND_MAX;
    }

    int max_iters = 100;
    int iter = 0;
    double tolerance = 1e-4;

    int num_threads = omp_get_max_threads();

    struct Cluster** local_clusters = malloc(num_threads * sizeof(struct Cluster*));
    for (int t = 0; t < num_threads; t++) {
        local_clusters[t] = malloc(K * sizeof(struct Cluster));
        for (int k = 0; k < K; k++) {
            local_clusters[t][k].points = malloc((n / num_threads + 1) * sizeof(struct Point));
            local_clusters[t][k].size = 0;
            local_clusters[t][k].capacity = n / num_threads + 1;
        }
    }

    while (iter < max_iters) {
        for (int i = 0; i < K; i++) {
            clusters[i].size = 0;
        }

        for (int t = 0; t < num_threads; t++)
            for (int k = 0; k < K; k++)
                local_clusters[t][k].size = 0;

#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            struct Cluster* my_clusters = local_clusters[tid];

#pragma omp for
            for (int i = 0; i < n; ++i) {
                struct Point point = points[i];
                int closest_index = 0;
                double min_distance = euclidean_distance(&point, &centroids[0]);

                for (int j = 1; j < K; ++j) {
                    double distance = euclidean_distance(&point, &centroids[j]);
                    if (distance < min_distance) {
                        min_distance = distance;
                        closest_index = j;
                    }
                }

                my_clusters[closest_index].points[my_clusters[closest_index].size++] = point;
            }
        }

        for (int t = 0; t < num_threads; t++) {
            for (int k = 0; k < K; k++) {
                int sz = local_clusters[t][k].size;
                memcpy(&clusters[k].points[clusters[k].size],
                       local_clusters[t][k].points,
                       sz * sizeof(struct Point));
                clusters[k].size += sz;
            }
        }

        double max_shift = 0.0;
#pragma omp parallel for reduction(max : max_shift)
        for (int i = 0; i < K; i++) {
            newCentroids[i] = calculate_centroid(
                clusters[i].points,
                clusters[i].size,
                centroids[i]);

            double shift = euclidean_distance(&centroids[i], &newCentroids[i]);
            if (shift > max_shift) max_shift = shift;
        }

        for (int i = 0; i < K; ++i) {
            centroids[i] = newCentroids[i];
        }

        if (max_shift < tolerance) break;

        iter++;
    }
    printf("Finished in %d iterations\n", iter);

    for (int t = 0; t < num_threads; t++) {
        for (int k = 0; k < K; k++)
            free(local_clusters[t][k].points);
        free(local_clusters[t]);
    }
    free(local_clusters);

    return clusters;
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