import numpy as np
import pytest
from initialize_kepler_cached import cached_cells


def test_resume_keeps_completed_initial_integrals(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path);sites=np.arange(10).reshape(5,2);cells=[np.array([[i,0],[i+1,0],[i,1]]) for i in range(5)];calls=[]
    def compute(i,site,cell):calls.append(i);return dict(index=i)
    assert cached_cells(sites,cells,2,compute,stop_after=2) is None
    assert calls==[0,1]
    rows=cached_cells(sites,cells,2,compute)
    assert calls==list(range(5));assert rows==[dict(index=i) for i in range(5)]
    cached_cells(sites,cells,2,lambda *args:pytest.fail('recomputed cached initial integral'))
    changed=sites.copy();changed[0,0]+=1
    with pytest.raises(ValueError,match='chunk mismatch'):cached_cells(changed,cells,2,compute)
