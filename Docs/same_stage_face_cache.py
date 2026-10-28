"""Explicit same-stage immutable face data; never shared between time stages."""
import hashlib,json
import numpy as np

class SameStageFaceCache:
    def __init__(self,stage_key):
        if not stage_key:raise ValueError('A bound state/geometry stage key is required')
        self.stage_key=stage_key
        self.fingerprint=None;self.cell_identity=None
        self.lines={};self.arcs={};self.line_hits=0;self.arc_hits=0

    def verify(self,cells,ids,positions,mass,factors,states,gen,priors):
        h=hashlib.sha256()
        for a in (ids,positions,mass,factors,gen):h.update(np.ascontiguousarray(a).view(np.uint8))
        for key in ('raw','R'):
            h.update(np.ascontiguousarray([s[key] for s in states]).view(np.uint8))
        for key in ('center','scale','projection','coefficient','standard_coefficient','residual_map'):
            h.update(np.ascontiguousarray([s['velocity_quadratic'][key] for s in states]).view(np.uint8))
        h.update(json.dumps(priors,sort_keys=True).encode())
        signature=h.hexdigest()
        if self.fingerprint is None:self.fingerprint=signature;self.cell_identity=id(cells)
        elif signature!=self.fingerprint or self.cell_identity!=id(cells):
            raise RuntimeError('A face cache cannot cross geometry or physical stages')

    def stats(self):
        return dict(stage_key=self.stage_key,immutable_fingerprint=self.fingerprint,
                    line_entries=len(self.lines),arc_entries=len(self.arcs),
                    line_hits=self.line_hits,arc_hits=self.arc_hits,
                    scope='Fixed geometry/velocity/prior only; each outer recomputes wn(z), roots and node physics')
