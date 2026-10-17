"""Create an experiment input or verify an identical input before resuming."""
import json
from pathlib import Path
from adaptive_coupled_mass import atomic


def install_input(path,payload):
    path=Path(path)
    if path.exists():
        previous=json.loads(path.read_text())
        if any(previous[key]!=payload[key] for key in ('config','t','state','mass','factors')):
            raise RuntimeError('continuation input changed')
        return False
    path.parent.mkdir(parents=True)
    atomic(path,payload)
    return True
