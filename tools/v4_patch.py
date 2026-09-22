#!/usr/bin/env python3
from pathlib import Path
import re, shutil

SRC=Path('MMS_BENCHMARKS_ABCD_TRUBA_V3')
ROOT=Path('MMS_BENCHMARKS_ABCD_TRUBA_V4')
if not SRC.exists():
    raise SystemExit('V3 generated root missing')
if ROOT.exists():
    shutil.rmtree(ROOT)
SRC.rename(ROOT)

for p in ROOT.rglob('*'):
    if p.is_file():
        try:
            s=p.read_text()
        except UnicodeDecodeError:
            continue
        s=s.replace('MMS_BENCHMARKS_ABCD_TRUBA_V3','MMS_BENCHMARKS_ABCD_TRUBA_V4')
        s=s.replace('V3_METHOD_CHANGES.md','V4_METHOD_CHANGES.md')
        p.write_text(s)

mg4=Path('tools/v4_mg_kernel.h').read_text()

solve_v4=r'''
static int solve(Sys*s,const char*solver,double tol,int maxit){
    int it=0;
    if(strcmp(solver,"SG_RBGS")==0){
        while(it<maxit && residual(s)>tol){ rbgs(s,2); it+=2; }
        return it;
    }
    while(it<maxit && residual(s)>tol){
        if(strcmp(solver,"MG2V")==0) mg4_cycle(s,2,2,1);
        else if(strcmp(solver,"MG2W")==0) mg4_cycle(s,2,2,2);
        else if(strcmp(solver,"MG3V")==0) mg4_cycle(s,3,2,1);
        else if(strcmp(solver,"RMT3H")==0) rmt_cycle_ref(s);
        else die("unknown solver");
        it+=1;
        if(!isfinite(residual(s))) break;
    }
    return it;
}
'''

for c in ROOT.glob('Benchmark_*/src/mms_solver.c'):
    s=c.read_text()

    # Benchmark D has an extra Qchem output parameter; patch the allocation
    # independently of the exact onevar_v2 signature.
    pat_grid=r'(static void onevar_v2\\([^\\n]*\\)\\{Sys s;)allocsys\\(&s,nx,ny\\);'
    if not re.search(pat_grid,s):
        raise SystemExit(f'V4 onevar grid anchor missing in {c}')
    s=re.sub(pat_grid,r'\\1allocsys(&s,nx+1,ny+1);',s,count=1)

    old="double cont=(BENCH_ID=='A')?projection_continuity(nx,ny,st,solver):0.0;"
    new="double cont=(BENCH_ID=='A')?projection_continuity(nx+1,ny+1,st,solver):0.0;"
    if old not in s: raise SystemExit(f'V4 projection grid anchor missing in {c}')
    s=s.replace(old,new,1)

    old='if(fieldout)write_real_field_v2(nx,ny,st,vars[0],solver,fieldout);'
    new='if(fieldout)write_real_field_v2(nx+1,ny+1,st,vars[0],solver,fieldout);'
    if old not in s: raise SystemExit(f'V4 field grid anchor missing in {c}')
    s=s.replace(old,new,1)

    marker='static int solve(Sys*s,const char*solver,double tol,int maxit)'
    pos=s.find(marker)
    if pos<0: raise SystemExit(f'V4 solve marker missing in {c}')
    s=s[:pos]+mg4+'\n'+s[pos:]

    pat=r'static int solve\(Sys\*s,const char\*solver,double tol,int maxit\).*?\nstatic void norms'
    if not re.search(pat,s,flags=re.S): raise SystemExit(f'V4 solve block missing in {c}')
    s=re.sub(pat,lambda _:solve_v4+'\nstatic void norms',s,flags=re.S)
    c.write_text(s)

for d in ROOT.glob('Benchmark_*'):
    bench=d.name.split('_')[1]
    slurm=d/f'benchmark_{bench}_array.slurm'
    slurm.write_text(f'''#!/bin/bash
#SBATCH --job-name=MMS_{bench}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --array=0-59%3
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err

set -euo pipefail
cd "${{SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}}"
GROUP_ID="${{SLURM_ARRAY_TASK_ID}}"
BASE_TASK=$((GROUP_ID * 5))
STATUS=0

echo "SLURM_ARRAY_JOB_ID=${{SLURM_ARRAY_JOB_ID:-NA}}"
echo "SLURM_ARRAY_TASK_ID=${{SLURM_ARRAY_TASK_ID}}"
echo "GROUP_ID=${{GROUP_ID}}"
echo "MANIFEST_TASKS=${{BASE_TASK}}-$((BASE_TASK + 4))"
echo "HOST=$(hostname)"

for OFFSET in 0 1 2 3 4; do
    TASK_ID=$((BASE_TASK + OFFSET))
    echo "Starting manifest task ${{TASK_ID}}"
    if ./run_case.sh "${{TASK_ID}}"; then
        echo "Finished manifest task ${{TASK_ID}}"
    else
        rc=$?
        echo "FAILED manifest task ${{TASK_ID}} rc=${{rc}}" >&2
        STATUS=1
    fi
done
if (( STATUS != 0 )); then exit 1; fi
echo "GROUP ${{GROUP_ID}} COMPLETE"
''')

    (d/'submit_all.sh').write_text(f'''#!/usr/bin/env bash
set -euo pipefail
mkdir -p logs
echo "Submitting Benchmark {bench}: 60 physical groups x 5 solvers = 300 executions"
echo "Maximum simultaneous physical groups = 3"
sbatch benchmark_{bench}_array.slurm
''')

    (d/'submit_task0.sh').write_text(f'''#!/usr/bin/env bash
set -euo pipefail
mkdir -p logs
sbatch --array=0 benchmark_{bench}_array.slurm
''')

    (d/'submit_index.sh').write_text(f'''#!/usr/bin/env bash
set -euo pipefail
group="${{1:?physical group index 0..59 required}}"
if (( group < 0 || group > 59 )); then
    echo "group must be 0..59" >&2
    exit 2
fi
mkdir -p logs
sbatch --array="${{group}}" benchmark_{bench}_array.slurm
''')
    for p in [slurm,d/'submit_all.sh',d/'submit_index.sh',d/'submit_task0.sh']:
        p.chmod(0o755)

note=r'''# MMS_BENCHMARKS_ABCD_TRUBA_V4

## Fair-comparison contract

V4 retains the validated V3 RMT structure and repairs the conventional
geometric-MG comparison.

All five methods solve the same fine-grid discrete equation, use the same zero
interior initial guess and the same residual stopping tolerance. The logical
campaign Nx x Ny values are interval/cell counts and the numerical system uses
Nx+1 by Ny+1 nodal points, so every campaign mesh is exactly nested under both
factor-two and factor-three coarsening.

MG2V/MG2W/MG3V use dimensionally consistent residual transfer: the integrated
fine residual is converted to residual density, tensor-product full weighting
is performed, then the coarse density is multiplied by coarse control-volume
area for the rediscretized integrated coarse equation. Coarse geometry is
exactly nested and prolongation is bilinear. There is no rollback, line search
or hidden SG-RBGS fallback.

RMT3H keeps the V3 production-structure adaptation: factor-three shifted
families, automatic independent x/y depth, one fine-index correction
workspace, coarsest-to-finest sawtooth, fine-defect CV restriction, direct tiny
coarsest solves, 16 noncoarsest postsmoothing sweeps and full correction.

Benchmark A is the direct pressure-solver comparison. B-D remain
operator-solver extensions.

## TRUBA submission

The 300-row manifest is preserved, but Slurm sees only 60 physical groups:

    #SBATCH --array=0-59%3

Each array task runs the five consecutive solver rows sequentially. This avoids
both MaxArraySize=100 and the previous 300-job/QOS submission problem while
retaining 300 total solver executions.
'''
(ROOT/'V4_METHOD_CHANGES.md').write_text(note)
for d in ROOT.glob('Benchmark_*'):
    (d/'docs'/'V4_FAIRNESS_AND_TRUBA_NOTE.md').write_text(note)

print('V4_PATCH_COMPLETE')
