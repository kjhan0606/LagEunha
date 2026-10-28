#!/bin/bash
# Box-side verification of Exam/Tests3D (no LagEunha build, no MPI job).
#   PY=/path/to/python-with-numpy-scipy-matplotlib ./selftest.sh
# 1. C reader/writer, 1D Evrard reference build
# 2. analytic solvers (Sedov xi0, Noh, KH linear theory vs Michalke)
# 3. every IC generator at small size + C-side totals check + byte round trip
# 4. every analyze.py --selftest on synthetic/analytic snapshots
# 5. every run.slurm: refuses a binary without the 3D marker; with a marker
#    and a stub mpirun it creates the run dir, IC and a fully filled params file
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
PY=${PY:-python3}
TMP=${TMP:-/tmp/lag3d_selftest}
rm -rf "$TMP"; mkdir -p "$TMP"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
fail=0
step() { echo; echo "=== $*"; }
res() { if [ "$1" = 0 ]; then echo "--> OK: $2"; else echo "--> FAILED: $2"; fail=1; fi; }

step "1. build C tools"
gcc -O2 -Wall -Wextra -std=c99 -o "$TMP/l3d_check" "$HERE/common/l3d_check.c" -lm; res $? "l3d_check.c"
gcc -O2 -Wall -o "$TMP/evrard1d_ref" "$HERE/Evrard3D/evrard1d_ref.c" -lm; res $? "evrard1d_ref.c"

step "2. analytic solvers"
$PY "$HERE/common/selftest_analytic.py"; res $? "analytic"

step "3. IC generators"
cd "$TMP"
ic() { local name=$1; shift; $PY "$@" -o "$TMP/$name.l3d" > "$TMP/$name.log" 2>&1; local r=$?;
       grep -E "np=|Ekin=|checks|FAILED|WARNING|hot|Noh:|Evrard:|KH3D:|wave:" "$TMP/$name.log";
       "$TMP/l3d_check" -w "$TMP/$name.c.l3d" "$TMP/$name.l3d" > "$TMP/$name.chk" || r=1
       tail -1 "$TMP/$name.chk"; cmp -s "$TMP/$name.l3d" "$TMP/$name.c.l3d" || { echo "C round trip differs"; r=1; }
       res $r "IC $name"; }
ic sw16      "$HERE/SoundWave3D/make_ic.py" --n 16
ic sw16boost "$HERE/SoundWave3D/make_ic.py" --n 16 --boost 1.0
ic sw16diag  "$HERE/SoundWave3D/make_ic.py" --n 16 --dir diag --mode displace
ic sed32     "$HERE/Sedov3D/make_ic.py" --n 32
ic sed16g    "$HERE/Sedov3D/make_ic.py" --n 16 --lattice glass
ic noh32     "$HERE/Noh3D/make_ic.py" --n 32
ic evr24     "$HERE/Evrard3D/make_ic.py" --n 24
ic evr16bg   "$HERE/Evrard3D/make_ic.py" --n 16 --background uniform
ic kh32      "$HERE/KH3D/make_ic.py" --n 32

step "4. analysis scripts on synthetic input"
for t in SoundWave3D Sedov3D Noh3D Evrard3D KH3D; do
	( cd "$HERE/$t" && $PY analyze.py --selftest --tmp "$TMP/$t" > "$TMP/$t.selftest.log" 2>&1 )
	r=$?; grep -E "SELFTEST|selftest:|exact:|1D ref|HydroCode|synthetic t" "$TMP/$t.selftest.log"; res $r "$t analyze.py --selftest"
done

step "5. run scripts (dry, stub mpirun, nothing submitted)"
mkdir -p "$TMP/fake/bin" "$TMP/fake/src"
ln -sfn "$(cd "$HERE/../.." && pwd)/Exam" "$TMP/fake/src/Exam"
printf '#!/bin/bash\necho fake\n' > "$TMP/fake/src/eunha2.exe"; chmod +x "$TMP/fake/src/eunha2.exe"
printf '#!/bin/bash\necho "stub mpirun $*"; exit 0\n' > "$TMP/fake/bin/mpirun"; chmod +x "$TMP/fake/bin/mpirun"
for t in SoundWave3D Sedov3D Noh3D Evrard3D KH3D; do
	case $t in SoundWave3D) n=16;; Sedov3D) n=16;; Noh3D) n=16;; Evrard3D) n=16;; KH3D) n=16;; esac
	LAGEUNHA_SRC="$TMP/fake/src" LAG3D_ROOT="$TMP/runs_nomarker" LAG3D_PY=$PY N=$n \
		bash "$HERE/$t/run.slurm" > "$TMP/$t.nomarker.log" 2>&1
	r1=$?
	cp "$TMP/fake/src/eunha2.exe" "$TMP/fake/src/eunha2.exe.bak"
	printf '#!/bin/bash\n# LAGEUNHA_3D_GFS_V1\necho fake\n' > "$TMP/fake/src/eunha2.exe"
	PATH="$TMP/fake/bin:$PATH" LAGEUNHA_SRC="$TMP/fake/src" LAG3D_ROOT="$TMP/runs" LAG3D_PY=$PY N=$n \
		bash "$HERE/$t/run.slurm" > "$TMP/$t.dry.log" 2>&1
	r2=$?
	mv "$TMP/fake/src/eunha2.exe.bak" "$TMP/fake/src/eunha2.exe"
	d=$(ls -d "$TMP/runs/${t}_"* 2>/dev/null | head -1)
	ok=0
	[ $r1 = 3 ] || { echo "no-marker exit $r1 (want 3)"; ok=1; }
	[ $r2 = 0 ] || { echo "dry exit $r2"; cat "$TMP/$t.dry.log"; ok=1; }
	[ -n "$d" ] && [ -s "$d/ic.l3d" ] && [ -s "$d/params3d.dat" ] && grep -q "EXIT:0" "$d/log" || { echo "run dir incomplete: $d"; ok=1; }
	[ -n "$d" ] && grep -q '@[A-Z_]*@' "$d/params3d.dat" && { echo "unfilled params"; ok=1; }
	[ -n "$d" ] && "$TMP/l3d_check" "$d/ic.l3d" > /dev/null || ok=1
	echo "$t: no-marker rc=$r1, dry rc=$r2, run dir $(basename "${d:-none}")"
	res $ok "$t run.slurm"
done

echo
if [ $fail = 0 ]; then echo "ALL SELFTESTS OK"; else echo "SOME SELFTESTS FAILED"; fi
exit $fail
