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


def _all_finite(*values):
    """True when every value is a finite float (excludes NaN and +/-inf).

    NaN/inf coordinates and box dimensions are not a physically meaningful
    input for these laws (SANITIZER.md 5.8's Precondition discipline: a
    checker's precondition must exclude genuinely undefined/ambiguous input
    classes), so callers use this to gate checkers whose invariant is only
    well-defined for finite inputs.
    """
    return all(math.isfinite(float(v)) for v in values)


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
    round trip for orthorhombic boxes, for any finite input box dimensions.

    Precondition excludes non-finite (NaN/inf) box dimensions: those are
    not a physically meaningful box and the round-trip law is not defined
    for them (SANITIZER.md 5.8's Precondition discipline).
    Tolerance is scaled by box length, same discipline as OM-PBC-001.
    """
    if not _all_finite(*set_dims, *got_dims):
        return
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

    Precondition excludes non-finite (NaN/inf) components: NaN is not
    equal to itself under IEEE 754, so an equality-based check of this
    exact identity is not meaningful for NaN inputs and is not a
    violation of the law (SANITIZER.md 5.8's Precondition discipline).
    No tolerance slack -- these are exact floating-point identities.
    """
    if not _all_finite(*v, *neg_neg_v, *v_plus_neg_v):
        return
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


# =====================================================================
# app/expandedensemblesampler.py, app/replicaexchangesampler.py --
# state-selection / exchange samplers
# =====================================================================

@_guard("expanded_ensemble_probability_normalization")
def check_probability_normalization(probability):
    """OM-SAMP-001: the log-sum-exp-normalized state probabilities in
    ExpandedEnsembleSampler.attemptStateChange must sum to 1.0.

    Tolerance scaled by the number of terms (standard floating-point
    summation error bound); empirically derived from a 200,000-trial
    sweep (see LAW_CANDIDATES.md OM-SAMP-001).
    """
    if not _all_finite(*probability):
        return
    n = max(1, len(probability))
    tol = 100 * eps64 * n
    trigger_if(not _within(sum(probability), 1.0, tol), "OM-SAMP-001")


@_guard("replica_exchange_criterion_consistency")
def check_exchange_criterion_consistency(exponent_uniform, ei_si_kt, ei_sj_kt, ej_sj_kt, ej_si_kt):
    """OM-SAMP-002: when ReplicaExchangeSampler.exchangeReplicas uses the
    uniform-kT Metropolis exponent formula, it must agree with the
    explicit-per-state-kT formula evaluated at a common temperature --
    two independently-coded expressions of the same exchange criterion.
    Arguments are already-reduced (dimensionless, divided by kT) energies,
    matching how the production code itself collapses Quantity/Quantity
    division to a plain float.

    Tolerance scaled by the largest reduced-energy term; empirically
    derived from a 200,000-trial sweep (see LAW_CANDIDATES.md OM-SAMP-002).
    """
    if not _all_finite(exponent_uniform, ei_si_kt, ei_sj_kt, ej_sj_kt, ej_si_kt):
        return
    exponent_explicit = (ei_si_kt - ej_si_kt) + (ej_sj_kt - ei_sj_kt)
    scale = max(abs(ei_si_kt), abs(ei_sj_kt), abs(ej_sj_kt), abs(ej_si_kt), 1e-12)
    tol = 1e6 * eps64 * scale
    trigger_if(not _within(exponent_uniform, exponent_explicit, tol), "OM-SAMP-002")


# =====================================================================
# amd.py -- aMD boost-energy formula vs. compiled force-scale factor
# =====================================================================

def _check_amd_effective_energy(energy, alpha, E, v_star, checker_id):
    """Shared body for OM-AMD-001/002: E - V* must equal
    alpha*(1 - alpha/(alpha+E-energy)) -- an exact algebraic identity
    relating the Python boost-energy formula to the literal force-scale
    sub-term (alpha/(alpha+E-energy)) embedded in the compiled engine's
    addComputePerDof expression string. Precondition: energy <= E (the
    boosted regime); for energy > E the compiled expression's own
    modify=step(E-energy) disables the boost on a different branch this
    law does not cover.

    Tolerance scaled by the larger energy-difference magnitude involved;
    empirically derived from a 200,000-trial sweep (see
    LAW_CANDIDATES.md OM-AMD-001).
    """
    if not _all_finite(energy, alpha, E, v_star) or energy > E:
        return
    denom = alpha + E - energy
    if denom == 0:
        return
    scale = alpha / denom
    lhs = E - v_star
    rhs = alpha * (1 - scale)
    tol = 1e6 * eps64 * max(1.0, abs(alpha), abs(E - energy))
    trigger_if(not _within(lhs, rhs, tol), checker_id)


@_guard("amd_total_energy_effective_energy")
def check_amd_total_effective_energy(energy, alpha, E, v_star):
    """OM-AMD-001: AMDIntegrator/DualAMDIntegrator's total-energy boost
    term. See _check_amd_effective_energy for the shared invariant.
    """
    _check_amd_effective_energy(energy, alpha, E, v_star, "OM-AMD-001")


@_guard("amd_group_energy_effective_energy")
def check_amd_group_effective_energy(energy, alpha, E, v_star):
    """OM-AMD-002: AMDForceGroupIntegrator/DualAMDIntegrator's
    force-group-energy boost term. See _check_amd_effective_energy for
    the shared invariant.
    """
    _check_amd_effective_energy(energy, alpha, E, v_star, "OM-AMD-002")


# =====================================================================
# app/metadynamics.py -- well-tempered metadynamics bias accumulation
# =====================================================================

@_guard("metadynamics_bias_accumulator_bound")
def check_metadynamics_bias_bound(total_bias_max, running_height_sum):
    """OM-META-001: the accumulated bias at any grid point must never
    exceed the sum of all Gaussian heights ever added, in the
    single-process case (Metadynamics.biasDir is None). Each Gaussian
    contribution to any single grid point is bounded by height*1.0 (the
    per-axis Gaussian kernel, and any outer product of per-axis kernels
    each in [0,1], peaks at 1.0), and _totalBias accumulates additively
    across calls, so the running maximum can never exceed the running
    sum of heights added.

    running_height_sum is an independent accumulator (SANITIZER.md 5.2):
    it is incremented once per _addGaussian call, alongside but not
    derived from the _totalBias array arithmetic itself -- comparing
    totalBias's own array against a quantity it was built from would be
    tautological; this compares it against a genuinely separate running
    total.

    Precondition (biasDir is None) excludes the multi-process case: when
    biases from other processes are loaded via _syncWithDisk, _totalBias
    incorporates those processes' own accumulated heights, which this
    process's running_height_sum does not track, so the bound would not
    hold in general for that case.
    """
    if not _all_finite(total_bias_max, running_height_sum):
        return
    tol = 1e6 * eps64 * max(1.0, abs(running_height_sum))
    trigger_if(total_bias_max - running_height_sum > tol, "OM-META-001")


# =====================================================================
# app/charmmparameterset.py -- NBFIX Lennard-Jones pair symmetry
# =====================================================================

@_guard("charmm_nbfix_reciprocity")
def check_nbfix_reciprocity(atom1_name, atom2_name, params_from_1, params_from_2):
    """OM-CHARMM-001: a CHARMM NBFIX term overrides the Lennard-Jones
    interaction between a specific pair of atom types. The Lennard-Jones
    pair potential is symmetric under particle exchange -- the
    interaction energy between atom type A and atom type B cannot depend
    on which one is labeled "first" -- so the stored per-atom-type NBFIX
    entries must be reciprocal: A's entry for B must equal B's entry for
    A (rmin, epsilon, rmin14, epsilon14).

    This is a genuine domain consequence, not internal bookkeeping:
    CharmmPsfFile.createSystem builds the CustomNonbondedForce
    acoef/bcoef tabulated-function matrix by looking up
    ``lj_type_list[i].nbfix[lj_type_list[j].name]`` independently for
    each ordered pair (i, j) (see charmmpsffile.py). A non-reciprocal
    NBFIX entry produces an asymmetric acoef/bcoef matrix, i.e. a
    simulated force where the energy of particle A pushing on B differs
    from B pushing on A -- an unphysical, non-Newtonian pairwise
    potential.

    params_from_1 / params_from_2 are the already-looked-up
    ``AtomType.nbfix`` dict entries (each a 4-tuple or None), read
    directly from production state -- this does not re-derive the NBFIX
    parsing logic, only compares two independently stored copies of what
    should be the same physical quantity (SANITIZER.md 5.2).

    Precondition: both entries must be present (an NBFIX line always
    calls add_nbfix on both atom types in the same statement, so by the
    time this observation point is reached after a successful NBFIX
    parse, both should already exist; a None here means the pair hasn't
    finished parsing yet, e.g. one atom type was undefined and silently
    skipped -- covered separately, not this law).
    """
    if params_from_1 is None or params_from_2 is None:
        return
    if not _all_finite(*params_from_1) or not _all_finite(*params_from_2):
        return
    trigger_if(tuple(params_from_1) != tuple(params_from_2), "OM-CHARMM-001")


@_guard("charmm_cmap_switch_range_bijection")
def check_cmap_switch_range_bijection(original_data, switched_data):
    """OM-CHARMM-002: _CmapGrid.switch_range() re-expresses a CMAP
    backbone-dihedral energy correction surface in a different angular
    coordinate convention (-180..180 degrees vs. 0..360 degrees), via a
    circular index shift. The physical correction-map surface itself --
    the set of energy values on the grid -- must not change; only the
    (phi, psi) angle labeling of each grid cell changes. switch_range
    must therefore be a bijection on grid values: every value present
    before must be present after, exactly once, with none created or
    dropped (SANITIZER.md 5.2's "stronger sanitizer" pattern: checking a
    semantic property -- the correction-map's physical content -- across
    a nontrivial transformation, not restating the transformation's own
    index arithmetic).

    This is a genuine domain consequence, not internal bookkeeping:
    CharmmPsfFile.createSystem feeds grid.switch_range().T directly into
    CMAPTorsionForce.addMap, so this transformed grid becomes the actual
    energy correction surface evaluated by the compiled engine at every
    simulation step for CMAP-corrected backbone dihedrals. A bijection
    violation (a duplicated or dropped energy value) would silently
    install a corrupted correction surface -- either double-counting an
    energy value at the expense of losing another the force field author
    intended, changing the shape of the backbone free-energy landscape.

    original_data / switched_data are the already-computed internal grid
    lists as production code stores them, read directly (not
    re-implementing the (i+mid)%res shift here).

    Precondition: none -- a grid of any resolution >= 1 is a valid CMAP
    input, and value-multiset preservation holds for every resolution
    (unlike a literal double-application involution check, which is only
    exact for even resolutions, since (res//2)*2 != res for odd res --
    that is a real mathematical fact about the truncating shift, not a
    defect, so it is deliberately not the law asserted here).
    """
    if not _all_finite(*original_data) or not _all_finite(*switched_data):
        return
    trigger_if(sorted(original_data) != sorted(switched_data), "OM-CHARMM-002")
