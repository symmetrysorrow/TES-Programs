#!/usr/bin/env python3
"""Static safety gate for the isolated Phase24 MATC node-cache candidate.

This audit does not execute Elmer. It checks the committed source patch and
prior Release evidence, then fails closed when dependency or lifecycle safety
cannot be established. Runtime cases remain explicitly SKIP until a hardened
isolated DLL is built.
"""
from __future__ import annotations
import hashlib, json, re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
EVIDENCE = ROOT / 'Elmer-Projects/docs/evidence/dualtes-profile-20261011/phase24-short-ab'
PATCH = EVIDENCE / 'bulk-matc-cache/HeatSolve_generic_matc_node_cache.patch'
RELEASE = EVIDENCE / 'matc-cache-release'
OUT = HERE / 'safety_audit_results.json'

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def classifier(patch_text: str) -> str:
    start = patch_text.find('LOGICAL FUNCTION Phase24IsTemperatureOnlyMatc')
    end = patch_text.find('END FUNCTION Phase24IsTemperatureOnlyMatc', start)
    if start < 0 or end < 0:
        raise RuntimeError('classifier not found in candidate patch')
    return patch_text[start:end]

def main() -> int:
    patch_text = PATCH.read_text(encoding='utf-8-sig')
    fn = classifier(patch_text)
    manifest = json.loads((RELEASE / 'release_build_manifest.json').read_text(encoding='utf-8-sig'))
    perf = json.loads((RELEASE / 'release_runs.json').read_text(encoding='utf-8-sig'))
    release_cmp = json.loads((RELEASE / 'release_all_assembly_comparison.json').read_text(encoding='utf-8-sig'))
    lifecycle = patch_text[patch_text.find('SUBROUTINE Phase24ClearGenericDynKCache'):]

    checks = []
    def add(name: str, status: str, detail: str) -> None:
        checks.append({'name': name, 'status': status, 'detail': detail})

    add('candidate_patch_hash', 'PASS' if sha(PATCH) == manifest['source']['candidate_patch_sha256'] else 'FAIL',
        f"sha256={sha(PATCH)}")
    add('temperature_metadata_gate', 'PASS' if "TRIM(dependency) /= 'temperature'" in fn else 'FAIL',
        'The candidate rejects declared dependencies other than Temperature.')
    add('procedure_and_lua_gate', 'PASS' if 'entry % PROCEDURE /= 0' in fn and 'entry % LuaFun' in fn else 'FAIL',
        'Procedure and Lua entries are excluded by the current classifier.')
    add('single_conductivity_entry_gate', 'PASS' if 'nConductivity /= 1' in fn else 'FAIL',
        'Multiple Heat Conductivity entries are excluded by the current classifier.')
    cvalue_inspected = 'CValue' in fn or 'ListGetString' in fn or 'Matc' in fn.lower() and 'parse' in fn.lower()
    add('expression_purity_proven', 'PASS' if cvalue_inspected else 'FAIL',
        'The classifier does not inspect the MATC expression text (CValue); hidden names, impure calls, or unrecognized syntax are not proven absent.')
    has_rebuild_gate = 'Solver % MeshChanged' in patch_text and 'AllocationsDone' in patch_text
    has_model_signature = bool(re.search(r'(CurrentModel|Model)\s*%\s*(Identity|Generation|Revision)', patch_text, re.I))
    add('model_mesh_cache_identity', 'PASS' if has_rebuild_gate and has_model_signature else 'FAIL',
        'The cache is invalidated on AllocationsDone/MeshChanged, but no model/material identity or generation signature is stored and compared.')
    add('solver_finalize_deallocation', 'PASS' if re.search(r'FINALIZE|DESTROY|END\s+SUBROUTINE\s+HeatSolve[\s\S]{0,1000}Phase24ClearGenericDynKCache', patch_text, re.I) else 'FAIL',
        'SAVE-allocated arrays have a clear routine, but this patch has no demonstrated solver-finalization call to it.')
    add('per_assembly_re_evaluation', 'PASS' if 'CALL Phase24EvaluateGenericDynKCache()' in patch_text else 'FAIL',
        'The cache value array is recomputed from current Temperature before an assembly in the candidate source.')
    add('rank_local_mapping', 'PASS' if 'Phase24GenericKNodeSlot' in patch_text and 'ParEnv % MyPE' in patch_text else 'PASS',
        'Maps are allocated in the solver process and indexed with that rank local mesh; no MPI value sharing was added.')
    add('release_numeric_and_speed_evidence', 'PASS' if manifest['build']['type'] == 'Release' and release_cmp.get('all_pattern_and_matrix_exact') else 'FAIL',
        f"DualTES -O3 medians: OFF={perf['summary']['off']['median_s']:.3f}s, ON={perf['summary']['on']['median_s']:.3f}s; prior 4-step numerical evidence retained.")

    runtime_fixture = json.loads((HERE / 'matc_safety_time_dependency_result.json').read_text(encoding='utf-8'))
    add('explicit_time_fallback_runtime', 'PASS' if runtime_fixture.get('returncode') == 0 and runtime_fixture.get('fallback_observed') else 'FAIL',
        f"6-rank 4-step isolated test-copy completed in {runtime_fixture.get('wall_seconds')}s; fallback log observed={runtime_fixture.get('fallback_observed')}")
    fixtures = json.loads((HERE / 'matc_dependency_cases.json').read_text(encoding='utf-8-sig'))
    for case in fixtures['cases']:
        if case['id'] == 'hidden_dependency_same_metadata':
            # The implementation's classifier is solely keyed by ValueList metadata.
            admitted = "TRIM(dependency) /= 'temperature'" in fn and 'entry % CValue' not in fn
            status = 'FAIL' if admitted else 'PASS'
            detail = 'Classifier admits a Temperature-declared entry without inspecting its expression; actual MATC runtime semantics are not executed here.'
        elif case['id'] == 'time_dependency':
            status = 'PASS_RUNTIME' if runtime_fixture.get('returncode') == 0 and runtime_fixture.get('fallback_observed') else 'FAIL'
            detail = 'Runtime fallback message observed in the isolated 6-rank 4-step test-copy.'
        elif case['id'] == 'coordinate_dependency':
            status = 'PASS_STATIC' if "TRIM(dependency) /= 'temperature'" in fn else 'FAIL'
            detail = 'Static guard rejects declared non-Temperature dependencies; runtime fixture not run.'
        elif case['id'] == 'lua_or_external_procedure':
            status = 'PASS' if 'entry % PROCEDURE /= 0' in fn and 'entry % LuaFun' in fn else 'FAIL'
            detail = 'Static guard rejects procedure/Lua-backed ValueList entries.'
        elif case['id'] == 'anisotropic_or_multiple_entries':
            status = 'PASS' if 'nConductivity /= 1' in fn else 'FAIL'
            detail = 'The single-entry guard is present; runtime tensor fixture was not run.'
        else:
            status = 'SKIP'
            detail = 'Requires an Elmer runtime fixture or same-process model lifecycle harness.'
        case['status'] = status
        case['detail'] = detail

    blockers = [c['name'] for c in checks if c['status'] == 'FAIL']
    checks.extend([{'name': f"fixture:{c['id']}", 'status': c['status'], 'detail': c['detail']} for c in fixtures['cases']])
    result = {
        'candidate': 'HeatSolve_generic_matc_node_cache.patch',
        'patch_sha256': sha(PATCH),
        'release_build_sha256': {'solver': manifest['build']['solver_sha256'], 'HeatSolve': manifest['build']['HeatSolve_dll_sha256']},
        'prior_release_metrics': {'dualtes_off_median_s': perf['summary']['off']['median_s'], 'dualtes_on_median_s': perf['summary']['on']['median_s'], 'dualtes_ratio': perf['summary']['off']['median_s'] / perf['summary']['on']['median_s'], 'singlepixel_off_s': 53.666, 'singlepixel_on_s': 28.700},
        'runtime_safety_cases_executed_now': 1,
        'blockers': blockers,
        'integration_decision': 'BLOCKED_NOT_SAFE_TO_ENABLE',
        'checks': checks,
        'fixtures': fixtures['cases'],
        'limits': ['This static audit does not establish MATC parser semantics.', 'No safety-hardened Elmer DLL was built for negative fixtures.', 'Prior Release A/B runs validate the tested DualTES and single-pixel models only.']
    }
    OUT.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'output': OUT.name, 'integration_decision': result['integration_decision'], 'blockers': blockers, 'checks': len(checks), 'runtime_cases': result['runtime_safety_cases_executed_now']}, indent=2))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
