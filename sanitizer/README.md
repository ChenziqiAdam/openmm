# Scientific sanitizer task: OpenMM (Reference platform)

This checkout is instrumented with **scientific sanitizers**: runtime checks that record an alarm when a
scientific or numerical law is violated. They never raise and never change a result. They are inactive unless the
environment variable `OPENMM_SANITIZE_LOG` names a file; an alarm appends one JSON line (`id`, `family`, `detail`).

## Your task

Write tests that trigger as many **distinct sanitizer families** as you can, through normal public OpenMM APIs
(the Python module `openmm`, Reference platform). Reading the sanitizer definitions and reasoning backwards from
their alarm conditions to inputs is part of the task.

- Definitions: `sanitizer/sanitizers.json` (precondition, invariant, observation point, alarm, rationale for each ID)
  and the instrumentation itself: `openmmapi/include/openmm/internal/Sanitizer.h` and the hook sites it names.
- Inputs must be valid uses of the API and satisfy each sanitizer's stated precondition. The Reference platform is the only
  one built.
- Primary metric: number of distinct triggered **families**. Secondary: distinct IDs, valid tests, cost.

## Environment

- Python with OpenMM installed: `~/Desktop/VSCode/openmm_venv/bin/python` (numpy and pytest available).
- The library and Python module are already built and frozen; you cannot rebuild them.
- Try a test and see which sanitizers fire:
  `OPENMM_SANITIZE_LOG=/tmp/my.jsonl ~/Desktop/VSCode/openmm_venv/bin/python -m pytest sanitizer/tests/test_x.py`
  then read `/tmp/my.jsonl`.

## Output

A test-only patch. Put pytest files in `sanitizer/tests/test_*.py` (data files next to them are allowed).
Tests must pass. Do **not**: edit production code or sanitizers, read or write `OPENMM_SANITIZE_*` variables or
the log from inside a test, call the sanitizer logger directly, use `subprocess`, `ctypes`, `os.environ`,
`importlib`, `eval`/`exec`, `PythonForce`/`CustomCPPForce`, or modify `sys.path`. Violations reject the submission.

## Scoring

The evaluator applies your patch to a fresh checkout of the frozen commit and runs each test twice (120 s each,
1800 s total). Only triggers present in both runs of a passing test count. Triggers are deduplicated to IDs and families.
