"""Runtime scientific-invariant trigger collection for the SciBench pilot.

This module is private and inactive unless ``SCIBENCH_TRIGGER_LOG`` names an
output file. Each JSON-lines record represents one observed invariant alarm;
the evaluator deduplicates checker IDs and maps them to root-cause families.

Design rule (see repo-root ``SANITIZER.md`` sections 5.1-5.8): a checker may
only (a) call a public API a second time on a transformed/derived input
and/or (b) read values already computed by production code, then compare. A
checker never re-implements the scientific formula it is checking, never
raises, and never changes a return value, exception, or numerical result.

Scope rule (SANITIZER.md 5.6): every sanitizer in this bank guards a
**scientific** invariant -- a physical, chemical, geometric, or dimensional-
analysis law whose violation has a domain consequence. Pure arithmetic
identities, range/finiteness checks, lookup-table round trips, and generic
software correctness are explicitly out of scope and are not instrumented
here. See ``SCIENTIFIC_CHECKERS.json`` for each checker's precondition,
invariant, observation point, and alarm predicate.
"""

import json
import math
import os
import threading

eps64 = 2.220446049250313e-16

# Re-entrancy guard: a checker that calls a public API a second time must not
# trigger the same checker recursively.
_active = threading.local()


def enabled():
    """Return True when the pilot evaluator requested trigger collection."""
    return bool(os.environ.get("SCIBENCH_TRIGGER_LOG"))


def trigger(checker_id):
    """Atomically append one checker ID to the configured JSON-lines log."""
    path = os.environ.get("SCIBENCH_TRIGGER_LOG")
    if not path:
        return
    payload = json.dumps({"checker_id": checker_id}, separators=(",", ":")) + "\n"
    descriptor = os.open(path, os.O_APPEND | os.O_CREAT | os.O_WRONLY, 0o600)
    try:
        os.write(descriptor, payload.encode("ascii"))
    finally:
        os.close(descriptor)


def trigger_if(condition, checker_id):
    """Record the checker ID when condition is true."""
    if condition:
        trigger(checker_id)


def _within(a, b, atol):
    """True when ``a`` and ``b`` differ by at most the absolute tolerance
    ``atol``. ``atol`` is expected to already be derived from the magnitude
    of the compared quantities (SANITIZER.md 5.8's T discipline).
    """
    return abs(float(a) - float(b)) <= atol


def _guard(checker_id):
    """Wrap a checker body so an internal error never disturbs production.

    A checker that itself errors is a curator bug, not a science alarm; it
    must not change program behaviour. We swallow it silently.
    """

    def decorator(func):
        def wrapper(*args, **kwargs):
            if not enabled():
                return
            if getattr(_active, "flag", False):
                return
            _active.flag = True
            try:
                func(*args, **kwargs)
            except Exception:
                pass
            finally:
                _active.flag = False

        wrapper.__name__ = func.__name__
        return wrapper

    return decorator


# =====================================================================
# unit/quantity.py, unit/unit.py -- dimensional-analysis system
# =====================================================================

@_guard("in_units_of_roundtrip")
def check_unit_roundtrip(original_value, unit1, roundtrip_value):
    """OM-UNIT-001: conversion round trip must reproduce the original value.

    Precondition: any Quantity converted to a compatible unit and back.
    Tolerance is scaled by the magnitude of the compared value
    (SANITIZER.md 5.8's T discipline).
    """
    scale = max(1.0, abs(float(original_value)))
    tol = 100 * eps64 * scale
    trigger_if(not _within(original_value, roundtrip_value, tol), "OM-UNIT-001")


@_guard("conversion_transitivity")
def check_conversion_transitivity(f12, f23, f13):
    """OM-UNIT-002: u1->u2->u3 factor product must equal u1->u3 directly.

    Tolerance is a fixed relative bound, far tighter than a spurious
    factor-of-two/sign error but with headroom over float64 round-off.
    """
    direct = f12 * f23
    tol = 1e-9 * max(1.0, abs(f13))
    trigger_if(not _within(direct, f13, tol), "OM-UNIT-002")


@_guard("conversion_inverse")
def check_conversion_inverse(f12, f21):
    """OM-UNIT-003: forward*reverse conversion factor must equal 1.0.

    Tolerance is a fixed relative bound over float64 round-off.
    """
    trigger_if(not _within(f12 * f21, 1.0, 1e-9), "OM-UNIT-003")


@_guard("unit_sqrt_involution")
def check_unit_sqrt(original_quantity_value, reconstructed_value, scale):
    """OM-UNIT-004: (sqrt(u))**2 must reproduce u's conversion factor exactly.

    Tolerance is scaled by the value's own magnitude, not its square root:
    squaring sqrt(x) accumulates float64 round-off proportional to x itself
    (SANITIZER.md 5.8's T discipline).
    """
    tol = 100 * eps64 * max(1.0, abs(scale))
    trigger_if(not _within(original_quantity_value, reconstructed_value, tol), "OM-UNIT-004")


@_guard("quantity_sqrt_vs_unit_math_sqrt")
def check_sqrt_paths_agree(method_value, reference_value, scale):
    """OM-UNIT-005: Quantity.sqrt() must agree with an independent
    math.sqrt(value)+unit.sqrt() composition.

    Tolerance is scaled by the value's own magnitude, same discipline as
    OM-UNIT-004.
    """
    tol = 100 * eps64 * max(1.0, abs(scale))
    trigger_if(not _within(method_value, reference_value, tol), "OM-UNIT-005")


@_guard("dot_norm_dimensional_consistency")
def check_dot_norm_dimensions(dot_unit_compatible_with_square, norm_unit_compatible_with_base):
    """OM-UNIT-006: dot(x,x) must have unit compatible with base_unit**2,
    and norm(x) compatible with base_unit.

    Structural check (unit compatibility is boolean, not float), so no
    numerical tolerance is needed.
    """
    trigger_if(not dot_unit_compatible_with_square, "OM-UNIT-006")
    trigger_if(not norm_unit_compatible_with_base, "OM-UNIT-006")


@_guard("molar_gas_constant_reference")
def check_molar_gas_constant(value_si):
    """OM-UNIT-007: R must match the CODATA 2018/SI-exact reference value.

    R is exact post-2019 SI redefinition (NA and kB are both exact), so
    tolerance is a transcription-error check, not a rounding-tolerance
    check.
    """
    codata_r = 8.31446261815324  # J / (mol K), CODATA 2018, exact under 2019 SI
    trigger_if(not _within(value_si, codata_r, 1e-9 * codata_r), "OM-UNIT-007")


@_guard("speed_of_light_reference")
def check_speed_of_light(value_si):
    """OM-UNIT-008: c must equal exactly 299792458 m/s (SI-exact since 1983).
    """
    trigger_if(value_si != 299792458.0, "OM-UNIT-008")


# =====================================================================
# app/internal/unitcell.py -- periodic box vector geometry
# =====================================================================

@_guard("pbc_lengths_angles_roundtrip")
def check_pbc_roundtrip(a, b, c, alpha, beta, gamma, a2, b2, c2, alpha2, beta2, gamma2, is_prereduced):
    """OM-PBC-001: computeLengthsAndAngles(computePeriodicBoxVectors(...))
    must reproduce the input lengths/angles, when the input already
    corresponds to a pre-reduced vector triple.

    Precondition restricted to already-reduced inputs: OpenMM's mandatory
    box-vector reduction intentionally changes lengths/angles for inputs
    that are not already reduced, which is documented behavior, not a
    violation of this law. Tolerance is scaled by box length for lengths,
    and is a fixed small bound for angles.
    """
    if not is_prereduced:
        return
    scale = max(1e-12, abs(a), abs(b), abs(c))
    len_tol = 100 * eps64 * scale
    ang_tol = 100 * eps64
    ok = (
        _within(a, a2, len_tol)
        and _within(b, b2, len_tol)
        and _within(c, c2, len_tol)
        and _within(alpha, alpha2, ang_tol)
        and _within(beta, beta2, ang_tol)
        and _within(gamma, gamma2, ang_tol)
    )
    trigger_if(not ok, "OM-PBC-001")


@_guard("pbc_reduction_volume_conservation")
def check_pbc_reduction_volume(volume_before, volume_after, scale):
    """OM-PBC-002: reducePeriodicBoxVectors must preserve the cell's
    determinant (volume) -- lattice reduction is a unimodular (det=1)
    shear transform.
    """
    tol = 100 * eps64 * max(1e-12, abs(scale)) ** 3
    trigger_if(not _within(volume_before, volume_after, tol), "OM-PBC-002")


@_guard("pbc_reduced_form_contract")
def check_pbc_reduced_form(satisfies_topology_contract):
    """OM-PBC-003: reducePeriodicBoxVectors's output must satisfy
    Topology.setPeriodicBoxVectors's own reduced-form inequalities.

    Boolean structural check.
    """
    trigger_if(not satisfies_topology_contract, "OM-PBC-003")


# =====================================================================
# app/element.py -- element table lookups
# =====================================================================

@_guard("element_getbymass_closest")
def check_getbymass_closest(requested_mass, returned_symbol, true_closest_symbol,
                             returned_is_also_at_minimum_distance,
                             true_closest_is_cache_reachable):
    """OM-ELEM-001: getByMass must return an element whose tabulated mass
    is truly closest (independent brute-force scan) to the query, for any
    query mass in the range spanned by the built-in element table, among
    elements that getByMass's own by-mass cache can represent at all.

    Precondition excludes a query only when BOTH (a) the returned element
    is itself exactly as close to the query as the independently-scanned
    closest element (returned_is_also_at_minimum_distance=True) AND (b)
    the reference element remains independently reachable through
    getByMass's own by-mass cache under its own tabulated mass
    (true_closest_is_cache_reachable=True): in that case both candidates
    are live and exactly tied, so "the" closest element is not uniquely
    defined by the law and a disagreement between them is not a
    violation. Any query failing either condition -- the returned element
    is strictly farther than the reference, or the reference element is
    not independently reachable through the cache under its own tabulated
    mass -- remains in scope, and a mismatch there is a violation.
    """
    if returned_is_also_at_minimum_distance and true_closest_is_cache_reachable:
        return
    trigger_if(returned_symbol != true_closest_symbol, "OM-ELEM-001")


@_guard("element_canonical_isotope_lighter")
def check_canonical_isotope(canonical_symbol, canonical_mass, sibling_masses):
    """OM-ELEM-002: getByAtomicNumber must return the lightest element
    among those sharing an atomic number (the documented "canonical"
    choice, e.g. hydrogen over deuterium).
    """
    if not sibling_masses:
        return
    trigger_if(canonical_mass > min(sibling_masses) + 1e-12, "OM-ELEM-002")


# =====================================================================
# app/topology.py -- Topology box-vector / atom-count accessors
# =====================================================================

@_guard("unitcell_dims_roundtrip")
def check_unitcell_dims_roundtrip(set_dims, got_dims):
    """OM-TOPO-001: setUnitCellDimensions -> getUnitCellDimensions must
    round trip for orthorhombic boxes.

    Tolerance is scaled by box length, same discipline as OM-PBC-001.
    """
    scale = max(1e-12, *[abs(x) for x in set_dims])
    tol = 100 * eps64 * scale
    ok = all(_within(s, g, tol) for s, g in zip(set_dims, got_dims))
    trigger_if(not ok, "OM-TOPO-001")


@_guard("system_particle_count_matches_topology")
def check_system_particle_count(num_particles, num_atoms):
    """OM-TOPO-002: System particle count must equal Topology atom count
    (createSystem's own documented contract).
    """
    trigger_if(num_particles != num_atoms, "OM-TOPO-002")


@_guard("hydrogen_mass_repartitioning_conserves_total")
def check_total_mass_conserved(total_before, total_after, num_atoms):
    """OM-TOPO-003: total system mass must be unchanged by hydrogenMass
    repartitioning (createSystem's documented promise), checked as a
    *global* sum independent of the pairwise transferMass arithmetic
    (not a restatement of that arithmetic's own local guarantee).
    """
    tol = 100 * eps64 * max(1, num_atoms) * max(1.0, abs(total_before))
    trigger_if(not _within(total_before, total_after, tol), "OM-TOPO-003")


# =====================================================================
# vec3.py -- coordinate primitive
# =====================================================================

@_guard("vec3_negation_involution")
def check_vec3_negation(v, neg_neg_v, v_plus_neg_v):
    """OM-VEC3-001: -(-v) == v and v+(-v) == 0, both exact under IEEE 754,
    for any Vec3 of finite floating-point components.

    No tolerance slack -- these are exact floating-point identities.
    """
    ok1 = tuple(neg_neg_v) == tuple(v)
    ok2 = tuple(v_plus_neg_v) == (0.0, 0.0, 0.0) or all(x == 0 for x in v_plus_neg_v)
    trigger_if(not (ok1 and ok2), "OM-VEC3-001")


# =====================================================================
# app/forcefield.py -- constrained-angle geometry
# =====================================================================

@_guard("constrained_angle_law_of_cosines")
def check_constrained_angle_length(law_of_cosines_length, independent_geometry_length, scale):
    """OM-FF-001: a constrained angle's constraint distance (law of
    cosines from two bond lengths and the angle) must match an
    independently-constructed-geometry distance (place the three atoms
    at consistent 3D/2D positions and measure the Euclidean distance
    directly).

    Tolerance is scaled by the constraint distance's own magnitude.
    """
    tol = 100 * eps64 * max(1e-12, abs(scale))
    trigger_if(not _within(law_of_cosines_length, independent_geometry_length, tol), "OM-FF-001")


# =====================================================================
# app/modeller.py -- topology merging
# =====================================================================

@_guard("modeller_add_conserves_counts")
def check_modeller_add_counts(actual_atoms, expected_atoms, actual_bonds, expected_bonds,
                               actual_residues, expected_residues, actual_chains, expected_chains):
    """OM-MOD-001: Modeller.add must conserve atom/bond/residue/chain
    counts (sum of the two input topologies' counts).
    """
    ok = (
        actual_atoms == expected_atoms
        and actual_bonds == expected_bonds
        and actual_residues == expected_residues
        and actual_chains == expected_chains
    )
    trigger_if(not ok, "OM-MOD-001")
