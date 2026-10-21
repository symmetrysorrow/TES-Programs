"""Ensure the installed native package includes the complete opt-in physics path."""

import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]


def test_setup_packages_required_physics_components():
    module = ast.parse((ROOT / "setup.py").read_text(encoding="utf-8"))
    build_class = next(n for n in module.body if isinstance(n, ast.ClassDef) and n.name == "CMakeBuild")
    run = next(n for n in build_class.body if isinstance(n, ast.FunctionDef) and n.name == "run")
    names_assign = next(n for n in ast.walk(run)
                        if isinstance(n, ast.Assign) and any(
                            isinstance(target, ast.Name) and target.id == "names" for target in n.targets
                        ))
    assert isinstance(names_assign.value, ast.IfExp)
    expected = {"dump2event", "posi2pulse", "dump2event_physical_reconstruct", "dumpall_transport_export"}
    for names in (ast.literal_eval(names_assign.value.body), ast.literal_eval(names_assign.value.orelse)):
        stems = {name.removesuffix(".exe") for name in names}
        assert stems == expected

    copied_files = [node for node in ast.walk(run)
                    if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "copy2"]
    assert any(isinstance(call.args[0], ast.Name) and call.args[0].id == "driver"
               for call in copied_files)
    assert (ROOT / "tes_cpp/dump2event/validation/build_physical_hdf5_v2.py").is_file()
