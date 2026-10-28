"""Add an isolated raw-rate diagnostic switch after the verified v5 patches."""

import hashlib
import json
from pathlib import Path
import sys


HEADER = "gfs_pressure_neutral_projection_20261027.h"
EXPECTED = "61b0b41f8c5ba9bb6321f5ad790626abc932dd25075675f15c924c7e782ce70c"
OLD = """    if (!count) return -3;
    if (!(norm>0.)) {"""
NEW = """    if (!count) return -3;
    /* Diagnostic only: retain the converged Jacobi rate without the global
     * pressure-work projection. The default branch below is unchanged. */
    const char *mode=getenv("GFS_GPU_VOLW2_PROJECTION");
    if (mode && strcmp(mode,"0")==0) {
        *alpha_out=0.;
        return 0;
    }
    if (mode && strcmp(mode,"1")!=0) return -5;
    if (!(norm>0.)) {"""


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Expected isolated source root")
    root = Path(sys.argv[1]).resolve()
    header = root / "Exam" / HEADER
    if digest(header) != EXPECTED:
        raise ValueError("Expected exact verified v5 projection header")
    source = header.read_text()
    if source.count(OLD) != 1 or source.count("#include <math.h>") != 1:
        raise ValueError("Projection header anchors changed")
    source = source.replace("#include <math.h>",
                            "#include <math.h>\n#include <stdlib.h>\n#include <string.h>")
    header.write_text(source.replace(OLD, NEW))
    print(json.dumps({
        "parent_header_sha256": EXPECTED,
        "patched_header_sha256": digest(header),
        "mode": "GFS_GPU_VOLW2_PROJECTION=0 bypasses global projection; "
                "unset or 1 retains the verified v5 operation",
    }, indent=2))


if __name__ == "__main__":
    main()
