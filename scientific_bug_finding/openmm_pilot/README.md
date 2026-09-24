# OpenMM scientific-sanitizer pilot

Fork: `https://github.com/ChenziqiAdam/openmm`, branch
`scibench-scientific-checkers-pilot`, built on pinned upstream tag `8.6.0`
(commit `c6173db6e8edd705eb59172bd21e9ce69c572405`). The fork's `master`
HEAD (`3c9effc9`) is 17 commits past the `8.6.0` release; `8.6.0` was
confirmed (via `git merge-base`) to be an ancestor of `master`, so the pin
uses the stable tagged release rather than the untagged dev tip, matching
this project's convention for the other pilots.

## Repository qualification (SANITIZER.md Step 0)

**Domain coverage**: OpenMM is a molecular dynamics engine (physics/
molecular simulation). This is a new domain for the bank -- the existing
pilots cover chemistry (pymatgen, pyscf), biology (biopython), astronomy
(astropy), geophysics/seismology (obspy), and quantum computing (qutip);
none previously covered classical/molecular-mechanics simulation.

**Engineering gate**: OpenMM has a mature CMake+CTest/pytest CI
(`.github/workflows/`), is a widely-used, citation-heavy (Eastman et al.,
PLOS Comp. Biol. 2017) MD engine with a stable public Python API, and its
Python layer (the `unit` dimensional-analysis system, the `app` module's
force-field/topology machinery) contains substantial, non-trivial
scientific logic of its own -- not merely a thin wrapper over the compiled
core (SANITIZER.md 5.6).

## Instrumentation location: pure-Python API layer, not the C++/CUDA/OpenCL core

OpenMM's scientific core (integrators, force evaluation, constraint
solvers) is implemented in C++ (`openmmapi/`, `platforms/`), with a
SWIG-generated Python wrapper. Building the full C++ core from source is a
heavy CMake+SWIG+platform-backend build. Before committing to an
instrumentation strategy, two things were verified:

1. A prebuilt `openmm` 8.6.1 conda package (`openmm_check` environment)
   resolves `openmm.unit.__file__`/`openmm.app.__file__` into
   `site-packages`, which was confirmed file-structure-identical (diff
   `-rq`, excluding compiled `.so`/`.pyx` artifacts) to this pilot's pinned
   checkout's `wrappers/python/openmm/unit/` and `app/` trees. Swapping
   these two pure-Python directories in the conda env's `site-packages`
   into the pinned checkout's copies gives a fully working, testable
   OpenMM installation without any C++/CUDA rebuild -- confirmed by
   running the bundled test suite (`TestForceField.py`, `TestElement.py`,
   `TestUnitCell.py`, `TestModeller.py`, `TestAPIUnits.py`,
   `TestTopology.py`, `TestUnits.py`) against the swapped installation.
2. Every scientific quantity examined in this pilot (energies, forces,
   positions/velocities, unit conversions, box geometry, element/topology
   bookkeeping) is either **computed entirely in the pure-Python layer**
   (`unit/`, `app/element.py`, `app/internal/unitcell.py`,
   `app/topology.py`/`forcefield.py`/`modeller.py`) or **queryable from
   Python via `Context.getState()`** after a compiled-core call returns --
   so no invariant in this bank requires instrumenting the compiled core
   itself.

**Decision: instrument only the pure-Python `unit`/`app`/`vec3.py` layer**
(Option A from the task brief), matching every other pilot's precedent
(their hooks live at Python call boundaries even when underlying numerics
are compiled elsewhere, e.g. pymatgen's Cython-backed neighbor search).
The compiled C++ core is untouched, so this pilot requires no rebuild
step for evaluation -- an agent only needs `pip install -e` against the
pinned checkout, exactly like the other pilots.

**Environment note**: verified `openmm.unit.__file__` /
`openmm.app.__file__` resolve into the pinned checkout's own
`wrappers/python/openmm/{unit,app}/` tree before any instrumentation or
testing began, following the same discipline the pymatgen pilot's README
documents.

## Candidate scan and law-first discipline

Per SANITIZER.md 5.7.1, `LAW_CANDIDATES.md` was written as a static
document -- precondition, invariant, observation point, alarm, rationale --
**before** any checker code was written or any candidate function was run
to observe its behavior. The scan covered `unit/quantity.py`, `unit/unit.py`,
`unit/unit_math.py`, `unit/constants.py`, `app/internal/unitcell.py`,
`app/element.py`, `app/topology.py`, `app/forcefield.py`'s bonded-force
generator classes, and `app/modeller.py`.

During Step 4 tolerance/precondition derivation (SANITIZER.md 5.8's T/X/P/N
discipline), two law candidates were found to be wrong as originally stated
and were corrected *before* instrumenting, not after:

- **OM-PBC-001** (box-vector <-> lengths/angles round trip): the naive
  precondition "any valid periodic box shape" was falsified by a direct
  sweep -- a triclinic input far from OpenMM's mandatory reduced-form
  convention produced a genuine 0.4-radian round-trip discrepancy that
  turned out to be *intentional, documented* reduction behavior, not a
  defect. The precondition was narrowed to inputs that already correspond
  to a pre-reduced vector triple; re-swept under the corrected precondition
  and found bit-level agreement.
- **OM-PBC-004** (originally: reduction preserves lengths/angles): hand-
  tracing the reduction arithmetic showed this claim is mathematically
  false in general (only `a`'s length is trivially invariant; `b`'s length
  and all three angles can and do change under reduction). Retired before
  instrumentation; superseded by the correctly-scoped OM-PBC-002 (volume
  conservation) and OM-PBC-003 (reduced-form contract), which are the
  provably-invariant properties of the same reduction step. See
  `LAW_CANDIDATES.md`'s rejected-candidates section for the full worked
  example of catching a wrong law before shipping it.

One law candidate (**OM-TOPO-003**, mass conservation under hydrogen-mass
repartitioning) was revised during formulation to avoid a SANITIZER.md 5.2
tautology: the per-pair `transferMass` arithmetic conserves the pairwise
sum by algebraic construction, so checking that in isolation would be a
restatement of the adjacent code. The accepted checker instead verifies
the *global* total-mass sum across all particles, which is not locally
guaranteed and would catch a bug in which atoms are selected for transfer.

## Bank summary

**19 sanitizers, 14 root-cause families**, spanning:

| subsystem | sanitizer count |
|---|---|
| `unit/` dimensional-analysis system (conversion round trips, transitivity, sqrt, dimensional consistency, physical constants) | 8 (OM-UNIT-001..008) |
| `app/internal/unitcell.py` periodic-box geometry | 3 (OM-PBC-001..003) |
| `app/element.py` element table lookups | 2 (OM-ELEM-001..002) |
| `app/topology.py`/`forcefield.py` System/Topology consistency | 3 (OM-TOPO-001..003) |
| `vec3.py` coordinate primitive | 1 (OM-VEC3-001) |
| `app/forcefield.py` constrained-angle geometry | 1 (OM-FF-001) |
| `app/modeller.py` topology merge count conservation | 1 (OM-MOD-001) |

This is one below the 20-sanitizer baseline in SANITIZER.md Step 0. The
scan was extended twice beyond the original `unit`/`app` core (into
`forcefield.py`'s bonded-force generators, then into `modeller.py`) before
accepting a below-baseline count rather than relaxing the quality bar; see
`LAW_CANDIDATES.md`'s summary for the specific further-scan targets that
were considered but not completed (charmm/amber/gromacs parameter-file
round-trips, `addSolvent`'s charge-neutralization step) -- an honest
next-step list for extending this bank further, not a claim that no more
sanitizers exist in OpenMM.

## Genuine finding: `Element.getByMass` duplicate-mass cache collision

While empirically deriving OM-ELEM-001's tolerance (a brute-force sweep of
every adjacent tabulated-mass midpoint), 8 genuine non-tie mismatches were
found between `getByMass`'s returned element and the true closest element.
Root cause, confirmed by reading `element.py`'s source and reproducing
directly against the unmodified library: `_elements_by_mass` is a plain
`dict` keyed by raw mass value, and two pairs of elements in the built-in
table share an exact tabulated mass (berkelium/curium both `247` Da;
dubnium/lawrencium both `262` Da -- both rounded whole-number estimates
for unstable synthetic elements). The dict construction silently drops the
first same-mass entry, so `getByMass` can **never** return curium or
lawrencium, for any input, even when one of them is the (tied-or-unique)
closest match.

This is a genuine OpenMM library defect, not a checker-tolerance artifact.
Searched `openmm/openmm` GitHub issues for "getByMass", "element mass
duplicate", "curium berkelium", "dubnium lawrencium" -- no existing report
found. Drafted (not filed, per this project's practice) as an ordinary
user bug report at `issues/ISSUE_1_getByMass_duplicate_mass_table_entries.md`,
with no mention of sanitizers/checkers/this benchmark project, per
SANITIZER.md 8.1.

## Verification performed

- **Isolated checker sensitivity** (SANITIZER.md Section 8): all 19
  checker functions in `_scientific_checkers.py` were called directly with
  synthetic bad states (a wrong roundtrip value, a mismatched transitivity
  product, a false compatibility flag, etc.) and every one recorded its
  expected checker ID. 19/19 pass.
- **Manual sanity exercises**: normal-usage calls through every
  instrumented entry point (`Quantity.in_units_of`, `conversion_factor_to`,
  `Quantity.sqrt`/`unit_math.sqrt`, `unit_math.norm`, module-level physical
  constants, `computePeriodicBoxVectors`/`reducePeriodicBoxVectors`,
  `Element.getByMass`/`getByAtomicNumber`, `Topology.setUnitCellDimensions`,
  `Vec3.__neg__`) with checkers enabled produced zero triggers.
- **Regression against the bundled test suite** (checkers enabled,
  `SCIBENCH_TRIGGER_LOG` set): `TestForceField.py`, `TestElement.py`,
  `TestUnitCell.py`, `TestModeller.py`, `TestAPIUnits.py`,
  `TestTopology.py`, `TestUnits.py` -- 197 tests, two rounds.
  **Baseline** (checkers disabled): 197 passed, 0 failed.
  **Round 1** (checkers enabled): 197 passed, 0 failed, but **2 checker
  triggers on legitimate library usage** -- `OM-TOPO-001` (in
  `TestForceField.py::test_PeriodicBoxVectors`) and `OM-UNIT-005` (in
  `TestUnits.py`, from `u.sqrt(1.0*kilogram*calorie)`). Both root-caused
  as checker-side bugs, not OpenMM defects (see "Checker-side fixes"
  below), and fixed. **Round 2** (post-fix, checkers enabled): 197 passed,
  0 failed, 0 triggers -- clean.
- A larger scenario exercising `createSystem` with `HAngles` constraints,
  `hydrogenMass` repartitioning, `addSolvent` (exercising periodic box
  construction), and `Modeller.add` (topology merge) on a real PDB
  structure (`ala_ala_ala.pdb` -> solvated system, 978 particles, 989
  constraints) produced zero triggers and no exceptions, both before and
  after the fixes.

## Checker-side fixes found by regression testing

Two checker bugs were found and fixed via the round-1 regression run
above -- both are the T/X/P/N-discipline failure modes SANITIZER.md 5.8
specifically warns about, caught here by running the checkers against the
library's own real test suite rather than only hand-written probes:

- **OM-UNIT-005** (unit-mismatch comparison): the original comparison
  squared the `sqrt` result and compared it against the input quantity's
  raw `._value` field without converting units first. `u.sqrt(1.0*
  kilogram*calorie)` genuinely differs numerically from its value in
  `kilogram*joule` (the calorie->joule factor, 4.184) despite being the
  same physical quantity -- a checker-side unit-handling bug. Fixed by
  converting the original quantity into the squared result's unit
  (`value_in_unit`) before comparing.
- **OM-TOPO-001** (missing side-channel invalidation): the checker-only
  side channel recorded by `setUnitCellDimensions` was never cleared when
  a later, unrelated `setPeriodicBoxVectors` call changed the same
  Topology's box vectors -- `TestForceField.py`'s `setUp` calls
  `setUnitCellDimensions(2,2,2)`, then `test_PeriodicBoxVectors` later
  calls `setPeriodicBoxVectors` with different vectors on the same
  instance without ever calling `setUnitCellDimensions` again, so
  `getUnitCellDimensions` compared the new (correct) box against a stale
  recording. Fixed by invalidating the side channel inside
  `setPeriodicBoxVectors`.

Both fixes were re-verified: the exact reproducing inputs are now silent,
a synthetic bad-state injection still fires (proving the checker did not
just go silent), and the full 197-test regression suite is clean across
round 2.

## Triggerability rounds

Two rounds were run, integrated with the tolerance-derivation sweeps
documented in `LAW_CANDIDATES.md` and above (the OM-PBC-001
precondition-narrowing sweep, the OM-ELEM-001 116-boundary brute-force
sweep, and the two regression-suite passes). Round 1 produced one
precondition fix (OM-PBC-001), one candidate retirement (OM-PBC-004), one
genuine upstream finding (OM-ELEM-001's `getByMass` defect), and two
checker-side bug fixes (OM-UNIT-005, OM-TOPO-001, both found via the
round-1 regression pass and detailed above). Round 2 (the post-fix
regression re-run) was clean: 0 new checker-side false positives. This
satisfies SANITIZER.md's audit-depth guidance ("continue rounds until two
consecutive rounds surface no new checker-side false positive") for the
regression-suite-based round; a broader hand-crafted adversarial fuzzing
round (beyond what the bundled test suite exercises) was not additionally
run in this pass -- recorded as a next step below, not claimed as done.

## What remains (honest next steps)

- A broader hand-crafted adversarial triggerability round beyond the
  bundled test suite (targeted fuzzing of unit conversions across more
  unit pairs, more box-vector shapes near the reduced-form boundary, more
  force-field templates with constrained angles) for additional
  confidence beyond the two regression-suite rounds already completed.
- Extending the bank toward the 20+ baseline via the further-scan targets
  noted in `LAW_CANDIDATES.md`'s summary (file-format parameter-set
  round-trips in `charmmparameterset.py`/`amberprmtopfile.py`/
  `gromacstopfile.py`; a carefully-scoped `addSolvent` charge-neutrality
  check).
- A possible follow-up "trajectory-level" checker class (see
  `LAW_CANDIDATES.md`'s rejected-candidates list) for integrator energy
  conservation and SETTLE/SHAKE constraint satisfaction, which require
  running the compiled core through actual integration steps and are out
  of scope for this pass's lightweight per-call checker style.
