import argparse
import csv
import random
from pathlib import Path

parser = argparse.ArgumentParser(
    description="Generate reproducible jittered 2D WLS verification clouds."
)
parser.add_argument("--n", type=int, required=True,
                    help="Number of points along each side.")
parser.add_argument("--output", type=Path, required=True,
                    help="Output CSV path.")
parser.add_argument("--seed", type=int, default=2412)
parser.add_argument("--jitter", type=float, default=0.15)
args = parser.parse_args()

if args.n < 3:
    raise SystemExit("n must be at least 3.")
if not (0.0 <= args.jitter < 0.5):
    raise SystemExit("jitter must be in [0, 0.5).")

rng = random.Random(args.seed)
spacing = 1.0 / (args.n - 1)
points = []

for j in range(args.n):
    for i in range(args.n):
        x = i * spacing
        y = j * spacing

        is_boundary = (
            i == 0 or i == args.n - 1 or
            j == 0 or j == args.n - 1
        )

        if not is_boundary:
            x += rng.uniform(
                -args.jitter * spacing,
                args.jitter * spacing
            )
            y += rng.uniform(
                -args.jitter * spacing,
                args.jitter * spacing
            )

        points.append((x, y))

args.output.parent.mkdir(parents=True, exist_ok=True)

with args.output.open("w", newline="") as f:
    writer = csv.writer(f)
    for x, y in points:
        writer.writerow([f"{x:.12f}", f"{y:.12f}"])

print(f"Generated: {args.output}")
print(f"Points: {len(points)} ({args.n} x {args.n})")
print(f"Spacing h: {spacing:.12f}")
print(f"Search radius for R=3.2h: {3.2 * spacing:.12f}")
print(f"Jitter fraction: {args.jitter}")
print("Domain: [0, 1] x [0, 1]")
