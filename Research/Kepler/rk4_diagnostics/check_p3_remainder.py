#!/usr/bin/env python3
"""80-digit independent algebra check of a native INITIAL basis, no evolution."""
import argparse
import csv
import json
from pathlib import Path
import mpmath as mp


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    destination = Path(args.output)
    assert not destination.exists()
    mp.mp.dps = 80
    alpha = [(3, 0, 0), (0, 3, 0), (0, 0, 3), (2, 1, 0), (1, 2, 0), (0, 2, 1), (0, 1, 2), (1, 0, 2), (2, 0, 1), (1, 1, 1)]
    report = {"scope": "native initial velocity interpolation residual and first derivatives, against direct 80-digit interpolation; no gravity or time evolution", "precision_digits": 80, "levels": {}}
    for level in (0, 1):
        maximum, points = mp.mpf(0), []
        with (Path(args.directory) / f"level{level}_remainder.csv").open() as source:
            rows = list(csv.DictReader(source))
        for row in rows:
            f = lambda key: mp.mpf(float(row[key]))
            triangle = [(f("ax"), f("ay")), (f("bx"), f("by")), (f("cx"), f("cy"))]
            (ax, ay), (bx, by), (cx, cy) = triangle
            det = (bx-ax)*(cy-ay)-(by-ay)*(cx-ax)
            x, y, eps = f("X"), f("Y"), f("eps")
            def velocity(xx, yy):
                omega = (xx*xx+yy*yy+eps*eps)**(-mp.mpf(3)/4)
                return [-yy*omega, xx*omega]
            nodal = []
            for powers in alpha:
                nx = sum(powers[k]*triangle[k][0] for k in range(3))/3
                ny = sum(powers[k]*triangle[k][1] for k in range(3))/3
                nodal.append(velocity(nx, ny))
            def residual(xx, yy, component):
                l1=((cy-ay)*(xx-ax)-(cx-ax)*(yy-ay))/det
                l2=(-(by-ay)*(xx-ax)+(bx-ax)*(yy-ay))/det
                bary = [1-l1-l2, l1, l2]
                value = velocity(xx, yy)[component]
                for powers, v in zip(alpha, nodal):
                    shape = mp.mpf(1)
                    for k in range(3):
                        for j in range(powers[k]):
                            shape *= (3*bary[k]-j)/(powers[k]-j)
                    value -= shape*v[component]
                return value
            exact = [residual(x, y, k) for k in (0, 1)]
            exact += [mp.diff(lambda xx: residual(xx, y, 0), x), mp.diff(lambda yy: residual(x, yy, 0), y), mp.diff(lambda xx: residual(xx, y, 1), x), mp.diff(lambda yy: residual(x, yy, 1), y)]
            actual = [f(k) for k in ("Rx", "Ry", "Rxx", "Rxy", "Ryx", "Ryy")]
            length = max(mp.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2) for a in triangle for b in triangle)
            scale = max([abs(z) for z in exact[:2]] + [length*abs(z) for z in exact[2:]])
            error = max([abs(a-b)/scale for a,b in zip(actual[:2],exact[:2])] + [length*abs(a-b)/scale for a,b in zip(actual[2:],exact[2:])])
            maximum = max(maximum,error)
            points.append({"triangle":int(row["triangle"]),"relative_value_and_length_scaled_derivative_error":float(error)})
        report["levels"][str(level)]={"queries":len(points),"max_normalized_error":float(maximum),"points":points}
        assert maximum < mp.mpf("1e-9"), "stable native remainder disagrees with independent input interpolation"
    destination.write_text(json.dumps(report,indent=2)+"\n")
    print({k:v["max_normalized_error"] for k,v in report["levels"].items()})


if __name__ == "__main__":
    main()
