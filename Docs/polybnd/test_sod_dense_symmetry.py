import json
from compare_sod_dense_symmetry import locate_step


def test_reports_store_endpoint_times(tmp_path):
    for i,end in enumerate((.1,.3,.6),1):
        (tmp_path/f'step_{i}_report.json').write_text(json.dumps(dict(t=end,dt=(.1,.2,.3)[i-1])))
    assert locate_step(tmp_path,.05,3)==1
    assert locate_step(tmp_path,.1,3)==1
    assert locate_step(tmp_path,.2,3)==2
    assert locate_step(tmp_path,.4,3)==3
    assert locate_step(tmp_path,.6,3)==3
