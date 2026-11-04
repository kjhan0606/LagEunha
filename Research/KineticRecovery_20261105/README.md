# Unlimited-gradient kinetic recovery admissibility

See ../../Docs/UNIFIED_HYDRO_KINETIC_RECOVERY_2026-11-05.md. The original native
C face-weighted fit evaluates exact positive cell states on uniform 1D meshes.
A Mach8 velocity jump yields Prec/P=-1.962962715714 at two cells for both N16
and N32. Constant and exact-affine-energy controls pass. The independent
closed expression matches to roundoff. Input energy/momentum are unchanged.
No limiter, floor, cap, reset, or reservoir hides the negative recovery.

This retires the specified unlimited-LSQ recovery. It does not reject every
limited reconstruction, every FV method, or the net-mass ALE transport
identities. No Gdot/Qdot, new RHS, spectrum, cold replay or PDE follows.
Source, native fit, commands and output are retained with SHA256 provenance.
