# Exact conservative assembly using a source cell and its outgoing overlaps

Implementation option within the existing common-refinement contract,
not a different PDE, kinetic reservoir or pressure correction. Not yet
implemented in the adaptive driver; that driver currently integrates all
fragments directly.

For each transported source cell S_i and target Laguerre partition C_j,
retain the original physical field and define D_ij=S_i intersect C_j.
For any transported moment f (including m,p,L,K,physical U), the exact identity is

    integral_(D_ii) f = integral_(S_i) f - sum_(j!=i) integral_(D_ij) f.

Therefore compute the source whole-cell moment once from the actual native
current state and integrate only outgoing overlaps with j!=i. Assemble the
target by incoming/outgoing moment transfers. Target mass remains
m_i_source-out_i+in_i, not a forced m_i_source. The approved mass/weight law
still needs an actual geometric solve. No target moment may be overwritten
after integration to manufacture agreement.

For this small finite deformation, outgoing fragments are thin. Whole-cell
K is already known to~2e-11 U by the independent native controls. Integrating
the small transferred K directly can avoid subtraction of two separately
quadrature-approximated large bulk-cell K values. It also avoids recomputing
the source bulk volume at each geometric iteration. This is an efficiency
and cancellation improvement; radial topology events must still be guarded.

The physical field representation must remain positive and geometric:
original source field restricted by the target half-plane predicates. The
own-cell fragment is the positive complement of its outgoing overlaps.
Do NOT turn the algebraic subtraction into negative physical quadrature
weights, an independent unresolved energy state or a heat top-up. Source
mass matrices minus actual outgoing-overlap mass matrices describe the
positive complement measure; verify positivity with integration uncertainty.

Source angular-boundary uncertainty in a region certified to remain owned
by i does not create intercell transfer and is already part of the known
whole-source integral. It should not be counted as an unknown outgoing flux.
Similarly, a radial chart panel where i uniformly beats every other site
has no transport-topology event. Pruning those cases is a sign certificate,
not a small-face cutoff. Keep all remaining unknown overlap bounds.

The926330-point preserved native source quadrature packet can independently
reconstruct source moments by ignoring its provisional target labels and
evaluating current chi,v. Future production should use the cached actual
reference mass/moment operators rather than reading an old test packet.
