from __future__ import annotations
import difflib, hashlib, re
from pathlib import Path
ROOT = Path(__file__).resolve().parent
BASE = Path(r'D:\Github\TES-Programs')
OLD_DUAL = ROOT.parent / 'dualtes-full' / 'case.sif'
OLD_STEADY = BASE / 'Elmer-Projects/workspace/generated/cases/case_toml_2cc534365acf_steady.sif'
RUN_STEADY = 'dualtes_wsl_steady_20261013_newton'
RUN_TRANSIENT = 'dualtes_wsl_transient_20261013_newton'
def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()
def common_settings(text: str) -> str:
    text = re.sub(r'(Nonlinear System Convergence Tolerance\s*=\s*)\S+', r'\g<1>1e-6', text)
    text = re.sub(r'(Linear System Max Iterations\s*=\s*)\S+', r'\g<1>4000', text)
    text = re.sub(r'(Linear System Convergence Tolerance\s*=\s*)\S+', r'\g<1>1e-12', text)
    text = re.sub(r'(BoomerAMG Relax Type\s*=\s*)\S+', r'\g<1>8', text)
    for key in ('Phase24 Dynamic K Node Cache', 'Phase24 Dynamic K Node Cache Verify', 'Phase24 Wall Profiling'):
        text = re.sub(r'^\s*' + re.escape(key) + r'\s*=.*\n', '', text, flags=re.M)
    for key, value in (
        ('Phase24 Vector Assembly', 'Logical True'),
        ('Phase24 Dynamic K Node Cache', 'Logical False'),
        ('Phase24 Dynamic K Node Cache Verify', 'Logical False'),
        ('Phase24 Precond Matrix Reuse', 'Logical True'),
        ('Phase24 Prism Lumped Preconditioner', 'Logical True'),
        ('Phase24 Coupled Temperature Tolerance', 'Real 1e-7'),
        ('TES Inner Circuit Newton', 'Logical True'),
    ):
        pattern = r'^\s*' + re.escape(key) + r'\s*=.*$'
        line = f'  {key} = {value}'
        if re.search(pattern, text, flags=re.M): text = re.sub(pattern, line, text, count=1, flags=re.M)
        else: text = re.sub(r'(  Variable DOFs\s*=\s*1\s*\n)', r'\1' + line + '\n', text, count=1)
    text = re.sub(r'(Steady State Convergence Tolerance\s*=\s*)\S+', r'\g<1>1e-6', text)
    return text
steady = OLD_STEADY.read_text(encoding='utf-8').replace('case_toml_2cc534365acf_steady', RUN_STEADY)
steady = common_settings(steady)
steady = re.sub(r'(Output File\s*=\s*).+', r'\1' + RUN_STEADY + '.result', steady)
(ROOT / 'inputs').mkdir(parents=True, exist_ok=True)
(ROOT / 'inputs' / f'{RUN_STEADY}.sif').write_text(steady, encoding='utf-8', newline='\n')
dual = OLD_DUAL.read_text(encoding='utf-8').replace('case_elongated_absorber_dual_tes_6a01a20b0b9e_pulse', RUN_TRANSIENT).replace('case_toml_2cc534365acf_steady', RUN_STEADY)
dual = common_settings(dual)
for key, value in (
    ('TES Inner Circuit Predict From Temperature', 'Logical True'),
    ('TES Inner Circuit Step Commit', 'Logical True'),
):
    pattern = r'^\s*' + re.escape(key) + r'\s*=.*$'
    line = f'  {key} = {value}'
    if re.search(pattern, dual, flags=re.M): dual = re.sub(pattern, line, dual, count=1, flags=re.M)
    else: dual = re.sub(r'(Variable DOFs\s*=\s*1\s*\n)', r'\1' + line + '\n', dual, count=1)
if not re.search(r'Phase24 BDF2 Predictor\s*=\s*Logical\s+True', dual): dual = re.sub(r'(BDF Order\s*=\s*2\s*\n)', r'\1  Phase24 BDF2 Predictor = Logical True\n', dual, count=1)
if not re.search(r'Phase24 BDF1 Step Start Times\(2\)', dual): dual = re.sub(r'(Phase24 BDF2 Predictor\s*=\s*Logical\s+True\s*\n)', r'\1  Phase24 BDF1 Step Start Times(2) = Real 0.02002 0.020020001\n', dual, count=1)
(ROOT / 'inputs' / f'{RUN_TRANSIENT}.sif').write_text(dual, encoding='utf-8', newline='\n')
for label, original, adapted, out in (('steady', OLD_STEADY, steady, ROOT/'inputs'/'steady_settings.diff'),('transient', OLD_DUAL, dual, ROOT/'inputs'/'transient_settings.diff')):
    diff=''.join(difflib.unified_diff(original.read_text(encoding='utf-8').splitlines(keepends=True),adapted.splitlines(keepends=True),fromfile=str(original),tofile=out.name))
    out.write_text(diff,encoding='utf-8',newline='\n')
    print(label,'source_sha256=',sha(original),'run_sha256=',hashlib.sha256(adapted.encode()).hexdigest(),'diff_lines=',len(diff.splitlines()))

