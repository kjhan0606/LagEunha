# Selected, qualified material profile. Source only from build scripts.
# Keeps the legacy Kepler/default compiler profile unchanged.
material_geometry_defines=(-DGFS_MATERIAL_COMMON_GEOMETRY)
material_host_defines=(
  -DGFS_MATERIAL_NATIVE_EXTERIOR
  -DGFS_MATERIAL_GRADIENT_PROJECT
  -DGFS_MATERIAL_PRESSURE_ACCELERATION
  -DGFS_MATERIAL_DIRECT_WALLS
  -DGFS_MATERIAL_POINT_GRAVITY
  -DGFS_MATERIAL_ORDERED_JUMP
  -DGFS_MATERIAL_FACE_NORMAL_LIMIT
  -DGFS_MATERIAL_THERMO_LIMIT_ONLY
  -DGFS_MATERIAL_DIRECT_FV
  -DGFS_MATERIAL_DIRECT_PAIR
  -DGFS_MATERIAL_DIRECT_PAIR_BATCH_BASIS
  -DGFS_MATERIAL_CENTROID_TRACE
  -DGFS_MATERIAL_CUDA_ACTIONS
  -DGFS_MATERIAL_BLAS_SCHUR
  -DGFS_MATERIAL_REUSE_RESPONSE
  -DGFS_MATERIAL_DIRECT_THREADED
  -DGFS_MATERIAL_REUSE_FACTOR
  -DGFS_MATERIAL_REUSE_GRAM
  -DGFS_MATERIAL_REUSE_RELEASE
  -DGFS_MATERIAL_CUDA_WEIGHT
  -DGFS_MATERIAL_CUDA_EXTERIOR
  -DGFS_MATERIAL_CUDA_RHS
  -DGFS_MATERIAL_CUDA_PROJECT
  -DGFS_MATERIAL_DIRECT_BASIS
  -DGFS_MATERIAL_ACTIVE_RELEASE
  -DGFS_MATERIAL_FACE_MEASURE_LS
  -DGFS_MATERIAL_COMMON_GEOMETRY
)
material_rhs_defines=(
  -DGFS_MATERIAL_GRADIENT_PROJECT
  -DGFS_MATERIAL_PRESSURE_ACCELERATION
  -DGFS_MATERIAL_DIRECT_WALLS
  -DGFS_MATERIAL_POINT_GRAVITY
  -DGFS_MATERIAL_FACE_NORMAL_LIMIT
  -DGFS_MATERIAL_DIRECT_FV
  -DGFS_MATERIAL_CUDA_ACTIONS_GATHER
)
# Device dispatch calls the same scalar projector as the CPU path.
material_project_defines=(
  -DGFS_MATERIAL_GRADIENT_PROJECT
  -DGFS_MATERIAL_PRESSURE_ACCELERATION
  -DGFS_MATERIAL_DIRECT_FV
  -DGFS_MATERIAL_DIRECT_WALLS
  -DGFS_MATERIAL_FACE_NORMAL_LIMIT
  -DGFS_MATERIAL_FACE_MEASURE_LS
)
# Same FP64 Jacobi-PCG and true-residual certificate, with iterations resident
# on the GPU. Qualified on saved KH/RT states and anchored native Cylinder
# evolution (200/400), with KH/RT/Sedov/Noh regressions; requires cooperative
# launch. Anchored solves retain the CPU original-edge true-residual check.
# Numerical history recycling remains opt-in, not part of this profile.
material_weight_defines=(-DGFS_MATERIAL_CUDA_PERSISTENT)
material_mkl_root=${MKLROOT:-/opt/ohpc/pub/intel/oneapi/mkl/latest}
material_host_includes=(-I"$material_mkl_root/include")
material_link_flags=(-L"$material_mkl_root/lib" -Wl,-rpath,"$material_mkl_root/lib" -lmkl_rt)
