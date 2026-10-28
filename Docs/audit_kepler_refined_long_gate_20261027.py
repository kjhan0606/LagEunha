"""Run the standard polar Kepler long-time gate for a chosen exact N."""

import argparse
from pathlib import Path
import sys

import audit_kepler_polar_long_gate_20261027 as standard_gate


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("run", type=Path)
    parser.add_argument("--expected-count", type=int, required=True)
    parser.add_argument("--time", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.expected_count <= 0 or args.time <= 0:
        raise ValueError("Expected a positive particle count and target time")

    # The shared gate is otherwise identical, but its baseline N is fixed at
    # 384^2. Override the module-level count so all snapshot and ID checks use
    # the refined 420^2 state without duplicating the scientific metrics.
    standard_gate.N = args.expected_count
    sys.argv = [str(Path(standard_gate.__file__)), str(args.initial),
                str(args.run), "--time", str(args.time),
                "--output", str(args.output)]
    standard_gate.main()


if __name__ == "__main__":
    main()
