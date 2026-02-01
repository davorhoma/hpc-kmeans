import random
import csv
import math
import argparse
import matplotlib.pyplot as plt


def generate_cluster(center_x, center_y, num_points, spread):
    points = []

    for _ in range(num_points):
        x = random.gauss(center_x, spread)
        y = random.gauss(center_y, spread)
        points.append((x, y))

    return points


def generate_dataset(
    clusters=3,
    points_per_cluster=100,
    spread=0.5,
    center_distance=5,
    seed=42
):
    random.seed(seed)

    all_points = []
    centers = []

    # Place cluster centers on a circle for good separation
    angle_step = 2 * math.pi / clusters

    for i in range(clusters):
        angle = i * angle_step
        cx = math.cos(angle) * center_distance
        cy = math.sin(angle) * center_distance

        centers.append((cx, cy))

        cluster_points = generate_cluster(
            cx, cy,
            points_per_cluster,
            spread
        )

        all_points.extend(cluster_points)

    return all_points, centers


def save_to_file(points, filename, csv_format=True):
    with open(filename, "w", newline="") as f:
        if csv_format:
            writer = csv.writer(f)
            for x, y in points:
                writer.writerow([x, y])
        else:
            for x, y in points:
                f.write(f"{x} {y}\n")


def plot_points(points, centers):
    xs, ys = zip(*points)
    cx, cy = zip(*centers)

    plt.scatter(xs, ys, s=10)
    plt.scatter(cx, cy, marker="x", s=100)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Generated K-means Dataset")
    plt.savefig("generator_clusters.png", dpi=300, bbox_inches="tight")
    print("Plot saved as generator_clusters.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--clusters", type=int, default=3)
    parser.add_argument("--points", type=int, default=100)
    parser.add_argument("--spread", type=float, default=0.5)
    parser.add_argument("--distance", type=float, default=5)
    parser.add_argument("--out", type=str, default="points.csv")
    parser.add_argument("--txt", action="store_true")
    parser.add_argument("--plot", action="store_true")

    args = parser.parse_args()

    points, centers = generate_dataset(
        clusters=args.clusters,
        points_per_cluster=args.points,
        spread=args.spread,
        center_distance=args.distance
    )

    save_to_file(points, args.out, not args.txt)

    print(f"Generated {len(points)} points")
    print("Cluster centers:", centers)

    if args.plot:
        plot_points(points, centers)
