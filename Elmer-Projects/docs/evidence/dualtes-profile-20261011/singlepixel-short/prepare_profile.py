from pathlib import Path
import json
import shutil
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from support.scenario_config import load_scenario, compile_project

scenario = load_scenario(ROOT / 'projects' / 'single_pixel.toml')
scenario['project'].setdefault('solver', {})['phase24_wall_profiling'] = True
model, mesh_name, steady_name, pulse_name = compile_project(scenario)
expected_mesh = 'mesh_hybrid_fullconf_h8'
if mesh_name != expected_mesh:
    raise RuntimeError(f'Unexpected mesh {mesh_name}; expected {expected_mesh}')
source_steady = 'case_toml_1b5ec032c418_steady'
case_name = 'case_singlepixel_profile_4step_20261011'
pulse = model['cases'][pulse_name]
pulse.update({
    'mesh': mesh_name,
    'restart_from': None,
    'restart_file_base': source_steady,
    'restart_time': 0.02,
    'preexisting_restart': True,
    'timesteps': [['18[us]', 1], ['1[us]', 2], ['1[ns]', 1]],
    'output_intervals': [1, 1, 1],
    'series_file': f'{case_name}_series.csv',
    'iteration_series_file': f'{case_name}_iterations.csv',
    'state_file': f'workspace/work/meshes/{mesh_name}/{case_name}.state',
    'output_result': False,
    'phase24_wall_profiling': True,
})
model['cases'] = {case_name: pulse}
out = Path(__file__).resolve().parent / 'profile_project.json'
out.write_text(json.dumps(model, indent=2) + '\n', encoding='utf-8')
mesh = ROOT / 'workspace' / 'work' / 'meshes' / mesh_name
src = mesh / f'{source_steady}.state'
dst = mesh / f'{case_name}.state'
if not src.is_file():
    raise FileNotFoundError(src)
shutil.copy2(src, dst)
print(f'case={case_name}')
print(f'mesh={mesh_name}')
print(f'four steps: 18 us, 1 us x2, 1 ns')
print(f'project={out}')
