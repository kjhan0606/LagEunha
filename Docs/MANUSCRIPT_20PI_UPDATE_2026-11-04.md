# Manuscript Kepler20pi update — 2026-11-04

The manuscript was updated and pushed to the author-specified Overleaf
project6ab774482cb870b10bb1a34f. Remote main is independently verified at
commit6cb0e1a471aa46aad65496a30ec23c5da5158dfa. The paper checkout is clean.
The title and bibliography are unchanged.

The abstract, test-method descriptions, Kepler section, discussion and
timestep appendix now record the exact20pi endpoint. The manuscript
distinguishes the RK4/globalKDK/two-levelKDK trajectory from a pure-integrator
test and preserves the separate HLLC standard-GFS and weighted-Q2 scopes.

Table1 presents the endpoint metrics and qualified published GIZMO velocity
comparison. Figure9 replaces the interim density figure with actualt0 and
20pi mass/area fields using linear .5..2 scaling. The caption and text
report density saturation and the measured tail. The original49185 disk
IDs,all90352 generators,core/buffer errors,full-t0 conservation and actual
production-path throughput are explicit. Same-time GIZMO density equivalence,
fullITS,multi-rank scaling and general eccentric-orbit accuracy are not claimed.
The remaining test list no longer lists20pi as unfinished.

Native snapshots were reread and both endpoint metrics recomputed. Timing
ratios were checked against original producer snapshot mtimes. The original
GADGET-2 eccentric-orbit RK4 statement was verified against the primary
paper. No citation keys or pre-existing labels were removed.

Local PDF is36pages with0overfull boxes,0undefined references/citations,
and0multiply-defined labels. Pages1,11,12,and13 were rendered and inspected.
The log retains4nonfatal float-placement warnings (the preceding build had3).
All updated figures and tables appear in the inspected PDF.

Authoritative source checkout:
/gpfs/kjhan/LagForce/overleaf_lageunha_publish_20261101_v1

Canonical local source:
/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs/lageunha.tex

Compiled PDF and full publication receipt:
/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/manuscript_20pi_update_20261104_v1/build/main.pdf
/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/manuscript_20pi_update_20261104_v1/publication_receipt.json

The Overleaf commit contains main.tex,the PDF/PNG density figure,and
data/kepler_20pi_metrics_20261104.json. Earlier figures and archives remain.
No simulation was launched,restarted,or modified for the manuscript update.
