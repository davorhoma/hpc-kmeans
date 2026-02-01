#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
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

    printf("input_file: %s\n", input_filename);
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

    printf("\n==== CLUSTERS ====\n");

    for (int i = 0; i < K; ++i) {
        printf("\nCluster %d (%d points):\n", i, clusters[i].size);
    }

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

    return 0;
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
        clusters[i].points = NULL;
        clusters[i].size = 0;
        clusters[i].capacity = 0;

        centroids[i].x = (double)rand() / RAND_MAX;
        centroids[i].y = (double)rand() / RAND_MAX;
    }

    int max_iters = 100;
    int iter = 0;
    double tolerance = 1e-4;

    while (iter < max_iters) {
        for (int i = 0; i < K; i++) {
            clusters[i].size = 0;
        }

        for (int i = 0; i < n; ++i) {
            struct Point point = points[i];
            int closestIndex = 0;
            double minDistance = euclidean_distance(&point, &centroids[0]);

            for (int j = 1; j < K; ++j) {
                double distance = euclidean_distance(&point, &centroids[j]);
                if (distance < minDistance) {
                    minDistance = distance;
                    closestIndex = j;
                }
            }

            add_point_to_cluster(&clusters[closestIndex], point);
        }

        double max_shift = 0.0;
        for (int i = 0; i < K; ++i) {
            struct Point newCentroid = calculate_centroid(clusters[i].points, clusters[i].size, centroids[i]);
            newCentroids[i] = newCentroid;

            double shift = euclidean_distance(&centroids[i], &newCentroid);
            if (shift > max_shift) max_shift = shift;

            centroids[i] = newCentroids[i];
        }

        if (max_shift < tolerance) break;

        iter++;
    }
    printf("Finished in %d iterations\n", iter);

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