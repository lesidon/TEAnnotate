#!/usr/bin/env python3
"""Check that a smoke test produced a usable software version report."""

import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("Version check requires PyYAML in the smoke-test Python environment.")

if len(sys.argv) != 2:
    sys.exit(f"Usage: {sys.argv[0]} SOFTWARE_VERSIONS_YML")

path = Path(sys.argv[1])
try:
    versions = yaml.safe_load(path.read_text(encoding="utf-8"))
except (OSError, yaml.YAMLError) as exc:
    sys.exit(f"{path}: {exc}")

if not isinstance(versions, dict) or not versions or any(
    not isinstance(tools, dict) or not tools for tools in versions.values()
):
    sys.exit(f"{path}: expected a nonempty mapping of steps to tool versions")

print(f"ok    {path} (parsed YAML, {len(versions)} version entries)")
