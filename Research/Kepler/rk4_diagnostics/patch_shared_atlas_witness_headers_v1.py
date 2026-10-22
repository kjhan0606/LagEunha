"""Add read-only hooks to a NEW frozen-run private header copy only."""
from pathlib import Path
import difflib
import sys

root=Path(__file__).resolve().parent
dest=Path(sys.argv[1]).resolve()
assert dest.parent==root and dest.name.startswith("shared_atlas_witness_")
changes={
    "material_shared_clip.h":[
        ("if(valid&&margin>best){best=margin;chosen=derivatives;out.dependent=dependent;}",
         "::atlas_witness::axis_trial(m,chart,dependent,active,derivatives,valid,margin);if(valid&&margin>best){best=margin;chosen=derivatives;out.dependent=dependent;}"),
        ("out.edge_windows+=roots.windows;need(roots.complete,roots.reason);events.insert(events.end(),roots.root.begin(),roots.root.end());",
         "out.edge_windows+=roots.windows;::atlas_witness::edge_events(index,boundary,roots);need(roots.complete,roots.reason);events.insert(events.end(),roots.root.begin(),roots.root.end());"),
        ("auto curve=boundary_curve(domain,edges,boundary,sf::range(Q(part)/16,Q(part+1)/16));for(int k=0;k<2;k++){lo[k]",
         "auto curve=boundary_curve(domain,edges,boundary,sf::range(Q(part)/16,Q(part+1)/16));::atlas_witness::radius_piece(m,boundary,part,curve);for(int k=0;k<2;k++){lo[k]"),
        ("need(r2>0,\"actual clipped reference region density-origin separation unproved\");return sf::lower(mf::root(S(r2)));",
         "::atlas_witness::radius_final(m,bound,lo,hi,r2);need(r2>0,\"actual clipped reference region density-origin separation unproved\");return sf::lower(mf::root(S(r2)));"),
    ],
    "material_shared_pullback.h":[
        ("need(step>0,\"real inverse iterate remains in source\");X-=step*dx;",
         "::atlas_witness::inverse_iteration(m,x,y,it,X,Y,fx,fy,dx,dy,step,lambda);need(step>0,\"real inverse iterate remains in source\");X-=step*dx;"),
        ("need(out.inside,\"real pullback inverse contraction box remains inside original source\");return out;",
         "::atlas_witness::inverse_final(m,x,y,X,Y,q,out,lambda,residual);need(out.inside,\"real pullback inverse contraction box remains inside original source\");return out;"),
    ],
}
diff=[]
for name,replacements in changes.items():
    path=dest/name
    original=path.read_text();modified=original
    for before,after in replacements:
        assert modified.count(before)==1,(name,before)
        modified=modified.replace(before,after)
    recovered=modified
    for before,after in reversed(replacements):recovered=recovered.replace(after,before)
    assert recovered==original
    path.write_text(modified)
    diff.extend(difflib.unified_diff(original.splitlines(True),modified.splitlines(True),fromfile="frozen_original/"+name,tofile="private_witness/"+name))
(dest/"READ_ONLY_OBSERVER_DIFF.patch").write_text("".join(diff))
