"""Unchanged physical stage with native straight-face and transport primitives."""
import hashlib
import importlib.util
import os
from pathlib import Path

_source=Path('/gpfs/kjhan/LagForce/Codex_full_disk_fast_modules_v1_20261019/full_disk_stage_memory.py')
_spec=importlib.util.spec_from_file_location('native_face_stage_implementation',_source)
_stage=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_stage)

from curved_face_nodes_native import build_geometry
from closed_root_transport_arrays import root_arrays
from closed_weight_root_native import solve_arrays
import native_straight_faces
import native_transport_root
import parallel_line_prepare

# Batch the unchanged per-point velocity formula across all straight-face
# samples owned by each cell, then feed the same native cut/fill implementation.
native_straight_faces.prepare=parallel_line_prepare.prepare

_stage._build_geometry=build_geometry
_stage.solve_from_transport=lambda geometry,transport,count:solve_arrays(root_arrays(geometry,transport),count)


def evaluate(*args,**kwargs):
    original_capture=_stage.capture
    journaling=os.environ.get('LAGFORCE_COUPLED_JOURNAL')=='1'
    if journaling:
        from full_moment_cell_journal import observe
        initial=args[0] if args else kwargs['initial']
        stage_key=args[6] if len(args)>6 else kwargs['stage_key']
        output=Path(args[8] if len(args)>8 else kwargs['output'])
        recorded=[]
        def capture(geometry,transport,values,states,gamma):
            captured=original_capture(geometry,transport,values,states,gamma)
            if not recorded:
                recorded.append(observe(geometry,transport,values,captured,initial['ids'],stage_key,initial['time'],output))
            return captured
        _stage.capture=capture
    try:result=_stage.evaluate(*args,**kwargs)
    finally:_stage.capture=original_capture
    # Keep the original passed audit report intact; publish added provenance
    # separately in the caller's complete-run report.
    result.native_provenance=dict(original_stage=str(_source),
        original_stage_sha256=hashlib.sha256(_source.read_bytes()).hexdigest(),
        straight_faces=native_straight_faces.provenance(),
        transport_root=native_transport_root.provenance(),
        batched_line_sample_preparation=parallel_line_prepare.provenance(),
        sparse_solver_changed=False,physical_closure_changed=False,
        observational_cell_journal=journaling)
    return result
