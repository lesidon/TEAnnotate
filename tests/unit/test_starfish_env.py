"""modules/local/starfish/main.nf: the shared Starfish environment is built once, under a
lock, and a ready environment is never deleted or replaced. The test runs the module's own
shell lines (between its starfish-env markers) with a fake solver."""
import re
import shutil
import subprocess
import time
from pathlib import Path

import pytest

MODULE = Path(__file__).resolve().parents[2] / 'modules' / 'local' / 'starfish' / 'main.nf'

pytestmark = pytest.mark.skipif(shutil.which('flock') is None or shutil.which('md5sum') is None,
                                reason='needs flock and md5sum (Linux)')

# fake micromamba: logs each build, takes a while, fails when FAIL_BUILD is set
SOLVER = '''#!/bin/bash
p=""; while [[ $# -gt 0 ]]; do [[ $1 == -p ]] && p=$2; shift; done
echo build >> "$BUILD_LOG"
sleep "${BUILD_SECONDS:-2}"
[[ -n "${FAIL_BUILD:-}" ]] && exit 1
mkdir -p "$p/bin"; printf '#!/bin/sh\\n' > "$p/bin/starfish"; chmod +x "$p/bin/starfish"
'''


def env_block(cache, lock_wait=None):
    block = re.search(r'# starfish-env-begin\n(.*?)\n\s*# starfish-env-end', MODULE.read_text(), re.S).group(1)
    for k, v in {"${params.starfish_env ?: ''}": '', '${cache_root}': str(cache),
                 '${env_yml}': 'env.yml', '${meta.id}': 'g'}.items():
        block = block.replace(k, v)
    if lock_wait is not None:
        block = block.replace('flock -w 10800 9', f'flock -w {lock_wait} 9')
    return block.replace('\\$', '$').replace('\\\\', '\\')


def start(tmp_path, name, cache, lock_wait=None, **env):
    """One task's environment step in its own directory, started in the background."""
    d = tmp_path / name
    (d / 'bin').mkdir(parents=True)
    (d / 'bin' / 'micromamba').write_text(SOLVER)
    (d / 'bin' / 'micromamba').chmod(0o755)
    (d / 'env.yml').write_text('name: starfish\n')
    script = 'set -ue\nOUT=out; mkdir -p out\n' + env_block(cache, lock_wait) + '\necho READY > out/status.txt\n'
    e = {'PATH': f"{d / 'bin'}:/usr/bin:/bin", 'BUILD_LOG': str(tmp_path / 'builds.log'), **env}
    return subprocess.Popen(['bash', '-c', script], cwd=d, env=e)


def status(tmp_path, name):
    return (tmp_path / name / 'out' / 'status.txt').read_text()


def envs(cache):
    return [p for p in cache.glob('starfish-*') if p.is_dir()]


def builds(tmp_path):
    return (tmp_path / 'builds.log').read_text().count('build')


def test_two_tasks_on_a_cold_cache_build_once(tmp_path):
    cache = tmp_path / 'cache'
    a, b = start(tmp_path, 'a', cache), start(tmp_path, 'b', cache)
    assert a.wait() == 0 and b.wait() == 0
    assert builds(tmp_path) == 1
    assert status(tmp_path, 'a') == status(tmp_path, 'b') == 'READY\n'
    (env,) = envs(cache)
    assert (env / '.teannotate_ok').exists() and (env / 'bin' / 'starfish').exists()


def test_a_failed_build_is_rebuilt_by_the_waiting_task(tmp_path):
    cache = tmp_path / 'cache'
    a = start(tmp_path, 'a', cache, FAIL_BUILD='1')
    time.sleep(0.5)                          # a holds the lock first
    b = start(tmp_path, 'b', cache, BUILD_SECONDS='0')
    assert a.wait() == 0 and b.wait() == 0
    assert status(tmp_path, 'a').startswith('skipped: building the Starfish environment failed')
    assert status(tmp_path, 'b') == 'READY\n' and builds(tmp_path) == 2
    (env,) = envs(cache)
    assert (env / '.teannotate_ok').exists()


def test_a_lock_timeout_skips_and_leaves_the_build_in_progress_alone(tmp_path):
    cache = tmp_path / 'cache'
    a = start(tmp_path, 'a', cache, BUILD_SECONDS='4')
    time.sleep(0.5)
    b = start(tmp_path, 'b', cache, lock_wait=1)
    assert b.wait() == 0 and status(tmp_path, 'b').startswith('skipped: timed out')
    assert a.wait() == 0 and status(tmp_path, 'a') == 'READY\n'
    (env,) = envs(cache)
    assert (env / '.teannotate_ok').exists() and builds(tmp_path) == 1
