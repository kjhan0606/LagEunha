"""Fixed typed wire contract for the native MPI-owned provider service."""
import ctypes as C
MAGIC=0x32564352534d464c
ABI=2
CONFIGURE,BOUNDS,DENSITY,VELOCITY,RELEASE,STOP=range(1,7)
FACE_STAGE,RESUME_STAGE=7,8
EXPORT_STATE,FACTOR_AUDIT=9,10
GEOMETRY,ADMISSIBILITY=11,12
DISPOSE_GEOMETRY=13
VERIFY_RESTART=14
COUPLED_METADATA=15
VERIFY_FORK,MIDPOINT_OBSERVE,COUPLED_JOURNAL=16,17,18
NATIVE_STAGE_REQUIRED=2
class Request(C.LittleEndianStructure):
    _pack_=1
    _fields_=[(k,C.c_uint64) for k in ('magic','epoch','run_token','request_id')]+[
        (k,C.c_uint32) for k in ('version','command','step','stage','outer','rank','nranks','ncells','ld_bytes','ld_mantissa')]+[
        ('physics_sha',C.c_char*64)]+[(k,C.c_char*512) for k in ('state_path','geometry_path','bounds_path','fits_path','output_dir','aux_path')]+[
        ('reserved',C.c_ubyte*888)]
class Response(C.LittleEndianStructure):
    _pack_=1
    _fields_=[(k,C.c_uint64) for k in ('magic','epoch','run_token','request_id')]+[
        (k,C.c_uint32) for k in ('version','status','ncells','ld_bytes','ld_mantissa','reserved')]+[
        (k,C.c_uint64) for k in ('payload_bytes','owned_start','owned_count','gravity_nodes','token')]+[
        ('output_path',C.c_char*512),('message',C.c_char*256),('padding',C.c_ubyte*160)]
assert C.sizeof(Request)==4096 and C.sizeof(Response)==1024
