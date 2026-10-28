"""Nested provider wall/process-CPU accounting without changing arithmetic.

Timer and wrapper overhead is included, not independently subtracted. No
warmup or calibration replay is performed. Use only in an isolated snapshot.
"""
from collections import defaultdict
from contextlib import contextmanager
from functools import wraps
from time import perf_counter_ns, process_time_ns

_stack = []
_rows = defaultdict(lambda: dict(calls=0, wall_ns=0, cpu_ns=0,
                               exclusive_wall_ns=0, exclusive_cpu_ns=0))


@contextmanager
def section(name):
    frame = [name, perf_counter_ns(), process_time_ns(), 0, 0]
    _stack.append(frame)
    try:
        yield
    finally:
        cpu = process_time_ns()-frame[2]
        wall = perf_counter_ns()-frame[1]
        assert _stack.pop() is frame
        row = _rows[name]
        row['calls'] += 1
        row['wall_ns'] += wall
        row['cpu_ns'] += cpu
        row['exclusive_wall_ns'] += wall-frame[3]
        row['exclusive_cpu_ns'] += cpu-frame[4]
        if _stack:
            _stack[-1][3] += wall
            _stack[-1][4] += cpu


def wrap(fn, name):
    @wraps(fn)
    def wrapped(*args, **kwargs):
        with section(name):
            return fn(*args, **kwargs)
    return wrapped


def report():
    assert not _stack
    return {name: dict(calls=row['calls'], **{key[:-3]+'_seconds': value/1e9
            for key, value in row.items() if key.endswith('_ns')})
            for name, row in sorted(_rows.items())}
