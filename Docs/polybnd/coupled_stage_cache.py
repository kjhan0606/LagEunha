"""Durable caches for completed numerical stages in a single trusted local run.

Pickle preserves ndarray dtypes and tuple structure. These files must only be
loaded from the run's own directory, never from untrusted inputs. Checksums
detect accidental corruption; they are not authentication.
"""
import hashlib
import json
import os
import pickle
from pathlib import Path


class StageCache:
    def __init__(self, directory, identity):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.identity = hashlib.sha256(json.dumps(identity, sort_keys=True,
            separators=(',', ':'), allow_nan=False).encode()).hexdigest()

    def __call__(self, name, compute):
        if name not in ('gravity_first', 'hydro_first', 'hydro_middle', 'gravity_last'):
            raise ValueError('unknown coupled stage')
        path = self.directory / (name + '.cache')
        if path.exists():
            with path.open('rb') as stream:
                header = json.loads(stream.readline())
                payload = stream.read()
            if header.get('version') != 1 or header.get('identity') != self.identity or header.get('stage') != name:
                raise ValueError('stage cache identity changed')
            if hashlib.sha256(payload).hexdigest() != header.get('sha256'):
                raise ValueError('stage cache checksum mismatch')
            return pickle.loads(payload)
        value = compute()
        payload = pickle.dumps(value, protocol=5)
        header = dict(version=1, identity=self.identity, stage=name,
            sha256=hashlib.sha256(payload).hexdigest())
        temporary = path.with_suffix('.tmp')
        with temporary.open('wb') as stream:
            stream.write(json.dumps(header, sort_keys=True).encode() + b'\n')
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        directory_fd = os.open(self.directory, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
        return value
