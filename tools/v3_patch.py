#!/usr/bin/env python3
from pathlib import Path
import re, shutil

SRC=Path('MMS_BENCHMARKS_ABCD_TRUBA_V2')
ROOT=Path('MMS_BENCHMARKS_ABCD_TRUBA_V3')
if not SRC.exists():
    raise SystemExit('V2 generated root missing')
if ROOT.exists():
    shutil.rmtree(ROOT)
SRC.rename(ROOT)

for p in ROOT.rglob('*'):
    if p.is_file():
        try:
            s=p.read_text()
        except UnicodeDecodeError:
            continue
        p.write_text(s.replace('MMS_BENCHMARKS_ABCD_TRUBA_V2',
                               'MMS_BENCHMARKS_ABCD_TRUBA_V3'))

kernel=Path('tools/v3_rmt_kernel.h').read_text()

solve_v3=r'''
static int solve(Sys*s,const char*solver,double tol,int maxit){
    int it=0;
    if(strcmp(solver,"SG_RBGS")==0){
        while(it<maxit&&residual(s)>tol){rbgs(s,2);it+=2;}
        return it;
    }
    while(it<maxit&&residual(s)>tol){
        if(strcmp(solver,"RMT3H")==0){
            rmt_cycle_ref(s);
            it+=1;
            if(!isfinite(residual(s)))break;
            continue;
        }
        double before=residual(s);
        size_t n=(size_t)s->nx*s->ny;
        double*save=malloc(n*sizeof(double));
        if(!save)die("save");
        memcpy(save,s->q,n*sizeof(double));
        if(strcmp(solver,"MG2V")==0){mg_cycle_ref(s,2,2,1);it+=46;}
        else if(strcmp(solver,"MG2W")==0){mg_cycle_ref(s,2,2,2);it+=86;}
        else if(strcmp(solver,"MG3V")==0){mg_cycle_ref(s,3,2,1);it+=46;}
        else {free(save);die("unknown solver");}
        double after=residual(s);
        if(!isfinite(after)||after>1.5*before){
            memcpy(s->q,save,n*sizeof(double));
            rbgs(s,4);
            it+=4;
        }
        free(save);
    }
    return it;
}
'''

for c in ROOT.glob('Benchmark_*/src/mms_solver.c'):
    s=c.read_text()
    pat=r'/\* V2 RMT production-flow correction:.*?\nstatic int solve\(Sys\*s,const char\*solver,double tol,int maxit\)'
    if not re.search(pat,s,flags=re.S):
        raise SystemExit(f'V2 RMT block not found in {c}')
    s=re.sub(pat,lambda _:kernel+'\nstatic int solve(Sys*s,const char*solver,double tol,int maxit)',s,flags=re.S)

    pat2=r'static int solve\(Sys\*s,const char\*solver,double tol,int maxit\).*?\nstatic void norms'
    if not re.search(pat2,s,flags=re.S):
        raise SystemExit(f'solve block not found in {c}')
    s=re.sub(pat2,lambda _:solve_v3+'\nstatic void norms',s,flags=re.S)
    c.write_text(s)

for d in ROOT.glob('Benchmark_*'):
    bench=d.name.split('_')[1]
    slurm=d/f'benchmark_{bench}_array.slurm'
    s=slurm.read_text()
    s=re.sub(r'^#SBATCH --array=.*\n','',s,flags=re.M)
    old='./run_case.sh "'+'$'+'{SLURM_ARRAY_TASK_ID}"'
    new='TASK_OFFSET="'+'$'+'{TASK_OFFSET:-0}"\n' \
        'TASK_ID=$((SLURM_ARRAY_TASK_ID + TASK_OFFSET))\n' \
        'echo "SLURM_ARRAY_JOB_ID='+'$'+'{SLURM_ARRAY_JOB_ID:-NA}"\n' \
        'echo "SLURM_ARRAY_TASK_ID='+'$'+'{SLURM_ARRAY_TASK_ID}"\n' \
        'echo "TASK_OFFSET='+'$'+'{TASK_OFFSET}"\n' \
        'echo "REAL_TASK_ID='+'$'+'{TASK_ID}"\n' \
        'echo "HOST=$(hostname)"\n' \
        './run_case.sh "'+'$'+'{TASK_ID}"'
    if old not in s:
        raise SystemExit(f'array dispatch line missing in {slurm}')
    slurm.write_text(s.replace(old,new))

    submit=f"""#!/usr/bin/env bash
set -euo pipefail
mkdir -p logs
echo "Submitting Benchmark {bench}: 300 cases as three 100-case arrays"
echo "TRUBA MaxArraySize=100"
echo "Maximum simultaneous cases = 3"

JOB1_RAW=$(sbatch --parsable --array=0-99%3 --export=ALL,TASK_OFFSET=0 benchmark_{bench}_array.slurm)
JOB1="${{JOB1_RAW%%;*}}"
JOB2_RAW=$(sbatch --parsable --dependency=afterany:${{JOB1}} --array=0-99%3 --export=ALL,TASK_OFFSET=100 benchmark_{bench}_array.slurm)
JOB2="${{JOB2_RAW%%;*}}"
JOB3_RAW=$(sbatch --parsable --dependency=afterany:${{JOB2}} --array=0-99%3 --export=ALL,TASK_OFFSET=200 benchmark_{bench}_array.slurm)
JOB3="${{JOB3_RAW%%;*}}"

cat > campaign_jobs.txt <<EOT
Benchmark {bench}
Batch 1: real tasks   0-99  -> job $JOB1
Batch 2: real tasks 100-199 -> job $JOB2
Batch 3: real tasks 200-299 -> job $JOB3
EOT
cat campaign_jobs.txt
"""
    (d/'submit_all.sh').write_text(submit)

    submit_idx=f"""#!/usr/bin/env bash
set -euo pipefail
idx="${{1:?real manifest index 0..299 required}}"
if (( idx < 0 || idx > 299 )); then
  echo "index must be 0..299" >&2
  exit 2
fi
offset=$(( (idx / 100) * 100 ))
local_id=$(( idx - offset ))
mkdir -p logs
sbatch --array="${{local_id}}" --export=ALL,TASK_OFFSET="${{offset}}" benchmark_{bench}_array.slurm
"""
    (d/'submit_index.sh').write_text(submit_idx)
    (d/'submit_task0.sh').write_text(
        f'#!/usr/bin/env bash\nset -euo pipefail\nmkdir -p logs\n'
        f'sbatch --array=0 --export=ALL,TASK_OFFSET=0 benchmark_{bench}_array.slurm\n'
    )
    for p in [slurm,d/'submit_all.sh',d/'submit_index.sh',d/'submit_task0.sh']:
        p.chmod(0o755)

note=r'''# MMS_BENCHMARKS_ABCD_TRUBA_V3

## RMT correction

V2 built an 18-way recursive allocation tree at every RMT level.  V3 removes
that tree and follows the RMT FINAL / FINALVOL2 computational structure: one
reusable fine-index correction workspace, automatic factor-three hierarchy,
one coarsest-to-finest sawtooth traversal, simultaneous shifted families,
fine-defect prefix reuse, direct tiny coarsest solves, 16 RBGS postsmoothing
sweeps per noncoarsest level, and full correction with omega=1.

There is no V2 residual-monotonic line search and no hidden SG-RBGS fallback in
the RMT branch.  Benchmark A is the direct pressure-solver comparison; B--D
remain operator-solver extensions.

## TRUBA MaxArraySize

The 300-row manifest is unchanged.  Because TRUBA MaxArraySize=100, submit_all.sh
creates three chained arrays with local indices 0--99 and offsets 0, 100, 200.
Each array retains %3 concurrency.  Slurm logs print TASK_OFFSET and REAL_TASK_ID.
'''
(ROOT/'V3_METHOD_CHANGES.md').write_text(note)
for d in ROOT.glob('Benchmark_*'):
    (d/'docs'/'V3_FIDELITY_AND_TRUBA_NOTE.md').write_text(note)

print('V3_PATCH_COMPLETE')
