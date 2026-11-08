import argparse
import random
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("n", type=int)
parser.add_argument("output")
parser.add_argument("--seed", type=int, default=2412)
parser.add_argument("--jitter", type=float, default=0.15)
args = parser.parse_args()

if args.n < 5:
    raise SystemExit("Use at least 5 points per axis.")

if not 0.0 <= args.jitter < 0.5:
    raise SystemExit("Jitter must be in [0, 0.5).")

n = args.n
h = 1.0 / (n - 1)
rng = random.Random(args.seed)

path = Path(args.output)
path.parent.mkdir(parents=True, exist_ok=True)

count = 0

with path.open("w") as file:
    for k in range(n):
        for j in range(n):
            for i in range(n):
                x = i * h
                y = j * h
                z = k * h

                if 0 < i < n - 1:
                    x += rng.uniform(-args.jitter * h, args.jitter * h)

                if 0 < j < n - 1:
                    y += rng.uniform(-args.jitter * h, args.jitter * h)

                if 0 < k < n - 1:
                    z += rng.uniform(-args.jitter * h, args.jitter * h)

                file.write(f"{x:.12f},{y:.12f},{z:.12f}\n")
                count += 1

print(f"Generated: {path}")
print(f"Grid: {n} x {n} x {n}")
print(f"Nodes: {count}")
print(f"Spacing: {h:.8e}")
print(f"Jitter: {args.jitter:.2f} H")
print(f"Seed: {args.seed}")
