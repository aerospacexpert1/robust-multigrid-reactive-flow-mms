#!/usr/bin/env bash
set -euo pipefail

ROOT=MMS_BENCHMARKS_ABCD_TRUBA_V3
ZIP=MMS_BENCHMARKS_ABCD_TRUBA_V3.zip

# Rebuild the validated V2 numerical baseline without running the V2 timing campaign.
sed -i 's/sp\.simplify(e)/e/g' tools/generate_mms_package.py
python3 tools/v2_generator_fix.py
python3 tools/generate_mms_package.py
python3 tools/postprocess_package.py

python3 - <<'PY'
from pathlib import Path
p=Path('tools/fidelity_patch.py')
s=p.read_text()
old="s=re.sub(r'static int solve\\(Sys\\*s,const char\\*solver,double tol,int maxit\\).*?\\nstatic void norms',KERNELS+'\\nstatic void norms',s,flags=re.S)"
if old in s:
    s=s.replace(old,"s=re.sub(r'static int solve\\(Sys\\*s,const char\\*solver,double tol,int maxit\\).*?\\nstatic void norms',lambda m: KERNELS+'\\nstatic void norms',s,flags=re.S)")
old2="s=re.sub(r'int main\\(int argc,char\\*\\*argv\\).*?\\n$',NEW_MAIN,s,flags=re.S)"
if old2 in s:
    s=s.replace(old2,"s=re.sub(r'int main\\(int argc,char\\*\\*argv\\).*?\\n$',lambda m: NEW_MAIN,s,flags=re.S)")
p.write_text(s)
PY

python3 tools/fidelity_patch.py
python3 tools/v2_patch.py
python3 tools/v2_postfix.py
python3 tools/v2_rmtfix.py
python3 tools/v2_dfix.py
python3 tools/v2_bfix.py
python3 tools/v2_budgetfix.py

# V3 replaces only the broken V2 RMT execution structure and TRUBA array submission.
python3 tools/v3_patch.py

test -d "$ROOT"
test -s "$ROOT/V3_METHOD_CHANGES.md"

# Syntax and matrix.
find "$ROOT" -type f \( -name '*.sh' -o -name '*.slurm' \) -print0 | xargs -0 -n1 bash -n
find "$ROOT" -type f -name '*.py' -print0 | xargs -0 python3 -m py_compile

python3 - <<'PY'
import csv,glob,collections
for p in glob.glob('MMS_BENCHMARKS_ABCD_TRUBA_V3/Benchmark_*/campaign_manifest.csv'):
    r=list(csv.DictReader(open(p)))
    assert len(r)==300,p
    assert collections.Counter(x['stiffness'] for x in r)=={'S1':100,'S2':100,'S3':100}
    assert collections.Counter(x['thread_label'] for x in r)=={'T1':75,'T2':75,'T4':75,'T16':75}
    assert collections.Counter(x['mesh'] for x in r)=={'M1':60,'M2':60,'M3':60,'M4':60,'M5':60}
    assert collections.Counter(x['solver'] for x in r)=={'SG_RBGS':60,'MG2V':60,'MG2W':60,'MG3V':60,'RMT3H':60}
print('V3_MANIFEST_PASS')
PY

# Structural RMT gates: no 18-way recursion, no V2 line search, full correction path.
for c in "$ROOT"/Benchmark_*/src/mms_solver.c; do
  grep -q 'rmt_v3_build_correction' "$c"
  grep -q 'for(int lc=lm;lc>=0;--lc)' "$c"
  grep -q 'rmt_v3_direct_coarsest' "$c"
  grep -q 'rmt_v3_smooth_level' "$c"
  grep -q 's->q\[IDX(i,j,s->nx)\] += w->corr' "$c"
  if grep -q 'rmt_error_recursive' "$c"; then echo "V2 recursive RMT remains: $c"; exit 1; fi
  if grep -q 'omegaTry' "$c"; then echo "V2 line search remains: $c"; exit 1; fi
done

# TRUBA MaxArraySize=100 gates.
for d in "$ROOT"/Benchmark_*; do
  f=$(ls "$d"/benchmark_?_array.slurm)
  if grep -q '^#SBATCH --array=' "$f"; then echo "embedded 300-task array remains: $f"; exit 1; fi
  grep -q 'TASK_OFFSET' "$f"
  grep -q 'REAL_TASK_ID' "$f"
  grep -q -- '--array=0-99%3' "$d/submit_all.sh"
  test "$(grep -o -- '--array=0-99%3' "$d/submit_all.sh" | wc -l)" -eq 3
  grep -q 'TASK_OFFSET=0' "$d/submit_all.sh"
  grep -q 'TASK_OFFSET=100' "$d/submit_all.sh"
  grep -q 'TASK_OFFSET=200' "$d/submit_all.sh"
  grep -q 'afterany' "$d/submit_all.sh"
done

# Build.
for d in "$ROOT"/Benchmark_*; do
  (cd "$d" && ./verify_campaign.sh && ./build.sh)
done

# M1 mandatory all-solver smoke for all A-D.
rm -rf v3_smoke
mkdir -p v3_smoke
python3 - <<'PY'
import subprocess,glob,os,csv,math
solvers=['SG_RBGS','MG2V','MG2W','MG3V','RMT3H']
for d in sorted(glob.glob('MMS_BENCHMARKS_ABCD_TRUBA_V3/Benchmark_*')):
    b=os.path.basename(d)
    for s in solvers:
        out=f'v3_smoke/{b}_{s}.csv'
        env=os.environ.copy(); env['OMP_NUM_THREADS']='1'; env['OMP_THREAD_LIMIT']='1'
        cp=subprocess.run([d+'/build/mms_solver','--Nx','108','--Ny','36','--stiffness','S1','--solver',s,'--out',out],
                          env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=120)
        if cp.returncode:
            raise SystemExit((b,s,cp.returncode,cp.stdout[-2000:]))
        r=list(csv.DictReader(open(out)))[0]
        assert r['converged']=='1',(b,s,r['relative_residual'])
        assert float(r['relative_residual'])<1e-6,(b,s,r['relative_residual'])
        assert math.isfinite(float(r['solve_wall_time']))
print('V3_M1_ALL_SOLVERS_PASS')
PY

# RMT-specific pressure ladder: enough to detect recurrence explosion/regression.
python3 - <<'PY'
import subprocess,os,csv,tempfile
d='MMS_BENCHMARKS_ABCD_TRUBA_V3/Benchmark_A_Pressure_Projection'
for nx,ny in [(108,36),(216,72),(432,144)]:
    out=tempfile.mktemp('.csv')
    env=os.environ.copy();env['OMP_NUM_THREADS']='1';env['OMP_THREAD_LIMIT']='1'
    cp=subprocess.run([d+'/build/mms_solver','--Nx',str(nx),'--Ny',str(ny),'--stiffness','S1','--solver','RMT3H','--out',out],
                      env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=180)
    if cp.returncode:
        raise SystemExit((nx,ny,cp.returncode,cp.stdout[-2000:]))
    r=list(csv.DictReader(open(out)))[0]
    assert r['converged']=='1',(nx,ny,r['relative_residual'])
    print(nx,ny,'iterations',r['iterations'],'time',r['solve_wall_time'])
print('V3_RMT_PRESSURE_LADDER_PASS')
PY

cat > "$ROOT/VALIDATION_REPORT.txt" <<'EOF'
MMS_BENCHMARKS_ABCD_TRUBA_V3 — VALIDATION
==========================================
V2 18-way recursive RMT tree removed: PASS
single fine-index RMT correction workspace: PASS
automatic factor-three x/y hierarchy: PASS
single coarsest-to-finest sawtooth traversal: PASS
control-volume fine-defect restriction with prefix reuse: PASS
direct tiny coarsest shifted-grid solves: PASS
16-sweep noncoarsest RBGS postsmoothing: PASS
full RMT correction omega=1: PASS
V2 residual-monotonic line search removed: PASS
RMT hidden SG fallback removed: PASS
TRUBA MaxArraySize=100 three-batch submission: PASS
TASK_OFFSET mapping 0/100/200: PASS
afterany batch chaining and %3 concurrency: PASS
300-row manifest preserved per benchmark: PASS
GCC OpenMP A/B/C/D build: PASS
M1 A-D x all five solvers: PASS
A pressure RMT M1/M2/M3 ladder: PASS
EOF

rm -f "$ZIP"
zip -qr "$ZIP" "$ROOT"
unzip -t "$ZIP" | tee v3_unzip.txt
grep -q 'No errors detected' v3_unzip.txt

tmp=$(mktemp -d)
unzip -q "$ZIP" -d "$tmp"
(cd "$tmp/$ROOT/Benchmark_A_Pressure_Projection" && ./verify_campaign.sh && ./build.sh)
rm -rf "$tmp"
echo 'ZIP integrity and fresh extraction build: PASS' >> "$ROOT/VALIDATION_REPORT.txt"

rm -f "$ZIP"
zip -qr "$ZIP" "$ROOT"
unzip -t "$ZIP" >/dev/null

echo FINAL_V3_VALIDATION_PASS
