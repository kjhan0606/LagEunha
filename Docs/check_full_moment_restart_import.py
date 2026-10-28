"""Test native checkpoint header rebind; no RHS or accepted physical step.

Use a real completed checkpoint. Numeric payloads and the old directory must
remain byte-identical. Model provenance verification is a separate prerequisite
for wiring this primitive into the evolution controller.
"""
import argparse
import ctypes as C
import hashlib
import json
from pathlib import Path


class Header(C.Structure):
    _fields_ = [(k, C.c_uint64) for k in (
        'magic', 'version', 'cells', 'epoch', 'run_token', 'state_elements',
        'ld_bytes', 'ld_mantissa')]+[('time', C.c_double), ('reserved', C.c_uint64*7)]


class State(C.Structure):
    _fields_ = [('header', Header)]+[(k, C.c_void_p) for k in (
        'ids', 'mass', 'state', 'factors', 'weight_rate_seed', 'Psi')]+[
        ('mapping', C.c_void_p*6), ('bytes', C.c_size_t*6)]


def digests(root):
    return {name: hashlib.sha256((root/(name+'.bin')).read_bytes()).hexdigest()
            for name in ('header', 'ids', 'mass', 'state', 'factors', 'Psi', 'weight_rate_seed')}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library', type=Path, required=True)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(exist_ok=False)
    lib = C.CDLL(str(args.library))
    lib.full_moment_state_import.argtypes = [C.c_char_p, C.c_char_p,
        C.POINTER(Header), C.c_uint64, C.POINTER(State)]
    lib.full_moment_state_import.restype = C.c_int
    lib.full_moment_state_close.argtypes = [C.POINTER(State)]
    lib.full_moment_state_close.restype = None
    original = digests(args.source)
    header = Header.from_buffer_copy((args.source/'header.bin').read_bytes())
    assert header.time > 0 and header.epoch and header.run_token
    new_token = header.run_token ^ 0x4000000000000000
    assert new_token and new_token != header.run_token
    result = State()
    destination = args.output/'rebound'
    assert lib.full_moment_state_import(str(args.source).encode(), str(destination).encode(),
        C.byref(header), new_token, C.byref(result)) == 0
    assert result.header.epoch == 1 and result.header.run_token == new_token
    expected = Header.from_buffer_copy(bytes(header))
    expected.epoch = 1
    expected.run_token = new_token
    assert bytes(result.header) == bytes(expected)
    assert (destination/'header.bin').read_bytes() == bytes(expected)
    lib.full_moment_state_close(C.byref(result))
    assert bytes(result) == bytes(C.sizeof(State))
    rebound = digests(destination)
    assert all(rebound[k] == value for k, value in original.items() if k != 'header')
    assert not (args.output/'rebound.pending').exists()

    failures = {}
    cases = [
        ('wrong_old_epoch', 'epoch', header.epoch+1, new_token),
        ('wrong_old_token', 'run_token', header.run_token+1, new_token),
        ('wrong_old_time', 'time', header.time+.001, new_token),
        ('incompatible_longdouble', 'ld_bytes', 8, new_token),
        ('same_run_token', None, None, header.run_token),
        ('zero_new_token', None, None, 0),
    ]
    for name, field, value, token in cases:
        target = args.output/name
        bad = Header.from_buffer_copy(bytes(header))
        if field is not None:
            setattr(bad, field, value)
        sentinel = State.from_buffer_copy(bytes([0x5a])*C.sizeof(State))
        before = bytes(sentinel)
        status = lib.full_moment_state_import(str(args.source).encode(), str(target).encode(),
            C.byref(bad), token, C.byref(sentinel))
        assert status and bytes(sentinel) == before
        assert not target.exists() and not Path(str(target)+'.pending').exists()
        failures[name] = True
    sentinel = State.from_buffer_copy(bytes([0x5a])*C.sizeof(State))
    before = bytes(sentinel)
    assert lib.full_moment_state_import(str(args.source).encode(), str(destination).encode(),
        C.byref(header), new_token+1, C.byref(sentinel)) != 0
    assert bytes(sentinel) == before and digests(destination) == rebound
    assert digests(args.source) == original
    report = dict(status='passed_native_checkpoint_header_import',
        source=str(args.source), time=header.time, cells=header.cells,
        payload_sha256={k: v for k, v in original.items() if k != 'header'},
        all_payload_bytes_exact=True, source_all_bytes_unchanged=True,
        only_epoch_and_token_changed=True, overwrite_rejected=True,
        invalid_import_no_output_change=failures,
        longdouble_bytes=header.ld_bytes, longdouble_mantissa=header.ld_mantissa,
        runtime_restart_connected=False, physical_providers=0, accepted_time_steps=0)
    (args.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
