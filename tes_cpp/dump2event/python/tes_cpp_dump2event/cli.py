"""Console entry point for the packaged native dump2event executable."""
from __future__ import annotations

import subprocess
import sys

from .converter import _executable


def main() -> int:
    """Run the native converter with the current command-line arguments."""
    return subprocess.run([_executable(), *sys.argv[1:]], check=False).returncode
