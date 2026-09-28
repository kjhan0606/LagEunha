"""I/O for the LagEunha 3D test suite (format "LAG3DV1").

One binary format is used for initial conditions *and* snapshots, so the
IC generators, the (future) 3D driver and the analysis scripts agree.
Little-endian, no padding. Layout:

  offset size  field
  0      8     magic          b"LAG3DV1\\0"
  8      4     int32 version  = 1
  12     4     int32 flags    bit0: pot[] present, bit1: vol[] meaningful
  16     8     int64 np
  24     8     double time
  32     8     double gamma
  40     48    double box[6]  xmin xmax ymin ymax zmin zmax
  88     12    int32 periodic[3]
  100    4     int32 reserved (0)
  104    24    double G, softening, reserved      (G = 0 for pure hydro)
  128    ...   SoA arrays, each np values:
               double x y z vx vy vz mass u rho vol   (u = SPECIFIC internal energy)
               [double pot]                           (if flags & 1)
               int64 id

The 2D drivers store ie = total internal energy of the cell (P V/(gamma-1)).
A 3D reader must set ie = mass*u. rho is the analytic/IC density (the driver
recomputes it from the Voronoi volume); vol is 0 in ICs unless flags&2.

The energy log that the analysis scripts expect from a 3D run is one line
per step of the form (key=value, order free, extra keys ignored):

  [E3D] step= 12 t= 0.0123 dt= 1e-3 Ekin= ... Eint= ... Epot= ... Etot= ...

This mirrors the 2D [RK4E] line. Epot is 0 when there is no gravity.
"""
import re
import struct
import numpy as np

MAGIC = b"LAG3DV1\0"
HEADER_SIZE = 128
FIELDS = ("x", "y", "z", "vx", "vy", "vz", "mass", "u", "rho", "vol")


def write_l3d(path, d, time=0.0, gamma=5.0 / 3.0, box=(0, 1, 0, 1, 0, 1),
              periodic=(1, 1, 1), G=0.0, softening=0.0, vol_valid=False):
    n = len(d["x"])
    flags = (1 if "pot" in d else 0) | (2 if vol_valid else 0)
    hdr = MAGIC + struct.pack("<iiqdd6d3ii3d", 1, flags, n, float(time),
                              float(gamma), *[float(b) for b in box],
                              *[int(p) for p in periodic], 0,
                              float(G), float(softening), 0.0)
    assert len(hdr) == HEADER_SIZE, len(hdr)
    with open(path, "wb") as f:
        f.write(hdr)
        for k in FIELDS:
            a = np.asarray(d.get(k, np.zeros(n)), dtype="<f8")
            assert a.shape == (n,), (k, a.shape)
            f.write(a.tobytes())
        if flags & 1:
            f.write(np.asarray(d["pot"], dtype="<f8").tobytes())
        ids = np.asarray(d.get("id", np.arange(n)), dtype="<i8")
        f.write(ids.tobytes())


def read_l3d(path):
    with open(path, "rb") as f:
        hdr = f.read(HEADER_SIZE)
        if hdr[:8] != MAGIC:
            raise ValueError("%s: not a LAG3DV1 file" % path)
        vals = struct.unpack("<iiqdd6d3ii3d", hdr[8:])
        version, flags, n, time, gamma = vals[:5]
        box = vals[5:11]
        periodic = vals[11:14]
        G, soft = vals[15], vals[16]
        d = {}
        for k in FIELDS:
            d[k] = np.frombuffer(f.read(8 * n), dtype="<f8").copy()
        if flags & 1:
            d["pot"] = np.frombuffer(f.read(8 * n), dtype="<f8").copy()
        d["id"] = np.frombuffer(f.read(8 * n), dtype="<i8").copy()
        if len(d["id"]) != n:
            raise ValueError("%s: truncated" % path)
    meta = dict(version=version, flags=flags, np=n, time=time, gamma=gamma,
                box=box, periodic=periodic, G=G, softening=soft)
    return d, meta


def density(d):
    """Snapshot density: rho if set, else mass/vol."""
    rho = d["rho"]
    if np.all(rho > 0):
        return rho
    vol = d["vol"]
    if np.all(vol > 0):
        return d["mass"] / vol
    raise ValueError("snapshot has neither rho nor vol")


def pressure(d, gamma):
    return (gamma - 1.0) * density(d) * d["u"]


def totals(d, G=0.0, pot=None):
    """Mass, momentum, kinetic, thermal energy (and potential if given)."""
    m = d["mass"]
    v2 = d["vx"] ** 2 + d["vy"] ** 2 + d["vz"] ** 2
    out = dict(M=m.sum(),
               Px=(m * d["vx"]).sum(), Py=(m * d["vy"]).sum(),
               Pz=(m * d["vz"]).sum(),
               Ekin=0.5 * (m * v2).sum(), Eint=(m * d["u"]).sum())
    if pot is not None:
        out["Epot"] = 0.5 * (m * pot).sum()
    return out


_KV = re.compile(r"([A-Za-z_][A-Za-z_0-9/|]*)=\s*([-+0-9.eEinfa]+)")


def parse_energy_log(path, tag="[E3D]"):
    """Parse '[E3D] key= value ...' lines into a dict of numpy arrays."""
    rows = []
    with open(path, errors="replace") as f:
        for line in f:
            if tag not in line:
                continue
            kv = {k: float(v) for k, v in _KV.findall(line.split(tag, 1)[1])}
            if "t" in kv:
                rows.append(kv)
    keys = sorted(set().union(*rows)) if rows else []
    return {k: np.array([r.get(k, np.nan) for r in rows]) for k in keys}


def periodic_delta(dx, L):
    return dx - L * np.round(dx / L)


def savefig(fig, path):
    fig.savefig(path, dpi=120, bbox_inches="tight")


def get_plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt
