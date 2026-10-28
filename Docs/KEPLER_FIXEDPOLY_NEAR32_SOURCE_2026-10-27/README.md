# Bounded Kepler near32 native source snapshot

This folder preserves the seven `Exam/` source/build files used for the
fixed-polygon near32 binary SHA-256
`3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e`.
The exact file hashes are in [MANIFEST.json](MANIFEST.json). Each copied
file was compared with the isolated native worktree on `syntax` after
copying. This is a source overlay, not a standalone Eunha source tree;
unchanged files come from the user's original source and the preserved
[native source overlay](../KEPLER_NATIVE_SOURCE_OVERLAY_2026-10-27/README.md).

The binary is preserved by the immutable A10 run directories, including
`/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/to0825_near32_v1/eunha2.exe`
and the ongoing `to65_near32_a10_v1` run. The bounded source adds an
opt-in fixed 1,024-gon wall, the W=0 energy-work correction, and exact
two-stage wall clipping. The near32 seed preclip is part of this
snapshot; it is still subject to later performance and geometry audits.

The [geometry-order report](../KEPLER_FIXEDPOLY_GEOMETRY_ORDER_2026-10-27.md)
records which early variants failed, the accepted-state replay, and the
first physical gate. This snapshot alone does not establish a stable
long-time disk or equivalence to the GIZMO paper's boundary treatment.
