"""Run replay.py unmodified, then once per named mutant.

The unmodified replay must exit 0 and print its PASS line. Every mutant must
make an assertion in replay.py fail: exit status 1 and a traceback ending in
AssertionError. A mutant that still passes means the corresponding check
would not notice that error.

Usage: python check_mutants.py
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import replay  # noqa: E402  (module-level code only defines names)

EXPECTED = 12


def run(*extra):
    return subprocess.run(
        [sys.executable, '-B', 'replay.py', *extra],
        cwd=ROOT, capture_output=True, text=True, timeout=300)


def main():
    if len(replay.MUTANTS) != EXPECTED:
        sys.exit(f'FAIL: replay.py lists {len(replay.MUTANTS)} mutants, expected {EXPECTED}')

    base = run()
    if base.returncode != 0 or 'PASS' not in base.stdout:
        sys.stderr.write(base.stdout + base.stderr)
        sys.exit(f'FAIL: unmodified replay exited {base.returncode}')
    print('baseline: exit 0,', base.stdout.strip().splitlines()[-1])

    bad = []
    for name in replay.MUTANTS:
        r = run('--mutant', name)
        lines = r.stderr.strip().splitlines()
        last = lines[-1] if lines else ''
        where = re.findall(r'replay\.py", line (\d+)', r.stderr)
        ok = r.returncode == 1 and last.startswith('AssertionError')
        print(f'mutant {name}: exit {r.returncode}, {last or "no error"}'
              + (f' at replay.py line {where[-1]}' if where else ''))
        if not ok:
            bad.append(name)
    if bad:
        sys.exit('FAIL: mutants not caught: ' + ', '.join(bad))
    print(f'PASS: baseline passes and all {EXPECTED} mutants are caught')


if __name__ == '__main__':
    main()
