import csv
import random
from pathlib import Path

# Reproducible point distribution
random.seed(2412)

# Domain and resolution
n = 21
spacing = 1.0 / (n - 1)
jitter_fraction = 0.15

output = Path("tests/data/wls_2d_verification.csv")
points = []

for j in range(n):
    for i in range(n):
        x = i * spacing
        y = j * spacing

        # Keep all four domain boundaries exactly fixed.
        is_boundary = (
            i == 0 or i == n - 1 or
            j == 0 or j == n - 1
        )

        if not is_boundary:
            x += random.uniform(
                -jitter_fraction * spacing,
                jitter_fraction * spacing
            )
            y += random.uniform(
                -jitter_fraction * spacing,
                jitter_fraction * spacing
            )

        points.append((x, y))

with output.open("w", newline="") as f:
    writer = csv.writer(f)
    for x, y in points:
        writer.writerow([f"{x:.12f}", f"{y:.12f}"])

print(f"Generated: {output}")
print(f"Number of points: {len(points)}")
print(f"Grid resolution: {n} x {n}")
print(f"Nominal spacing: {spacing:.6f}")
print(f"Maximum jitter per coordinate: {jitter_fraction * spacing:.6f}")
print(f"Domain: [0, 1] x [0, 1]")
