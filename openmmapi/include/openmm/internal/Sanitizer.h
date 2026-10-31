#ifndef OPENMM_SANITIZER_H_
#define OPENMM_SANITIZER_H_

/*
 * Scientific sanitizer support (observation only).
 *
 * Inactive unless the environment variable OPENMM_SANITIZE_LOG names a file.  An alarm appends one JSON line.
 * Checks read state that already exists at their observation point; they never modify it.
 */

#include "openmm/System.h"
#include "openmm/Vec3.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace OpenMM {
namespace sanitizer {

static const double EPS = 2.220446049250313e-16;

inline const char* logPath() {
    static const char* path = std::getenv("OPENMM_SANITIZE_LOG");
    return path;
}

inline bool enabled() {
    const char* p = logPath();
    return p != NULL && p[0] != '\0';
}

inline void alarm(const char* id, const char* family, const std::string& detail) {
    if (!enabled())
        return;
    FILE* f = std::fopen(logPath(), "a");
    if (f == NULL)
        return;
    const char* key = std::getenv("OPENMM_SANITIZE_KEY");
    std::fprintf(f, "{\"id\":\"%s\",\"family\":\"%s\",\"key\":\"%s\",\"detail\":\"%s\"}\n", id, family, key == NULL ? "" : key, detail.c_str());
    std::fclose(f);
}

inline std::string num(double a, double b = NAN) {
    char buf[96];
    if (std::isnan(b))
        std::snprintf(buf, sizeof(buf), "%.6g", a);
    else
        std::snprintf(buf, sizeof(buf), "%.6g vs %.6g", a, b);
    return buf;
}

inline double norm(const Vec3& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/**
 * Reaction-field shift: the shifted pair potential 1/r + krf r^2 - crf must vanish at the cutoff.
 */
inline void checkReactionField(double cutoff, double krf, double crf) {
    if (!enabled())
        return;
    // For a dielectric constant so large that 3*eps overflows the constants are not representable; there is no
    // potential to check.
    if (!std::isfinite(krf) || !std::isfinite(crf))
        return;
    double residual = 1.0/cutoff + krf*cutoff*cutoff - crf;
    if (std::fabs(residual) > 16*EPS/cutoff)
        alarm("RF-001", "rf_continuity", "potential at cutoff " + num(residual, 16*EPS/cutoff));
}

/**
 * A switching function takes values in [0,1] and does not increase with distance.
 */
inline void checkSwitch(double value, double deriv, double cutoff, double switching) {
    if (!enabled())
        return;
    double slack = 16*EPS;
    if (value < -slack || value > 1+slack || deriv > slack/(cutoff-switching))
        alarm("SW-001", "switching_bounds", "value " + num(value, deriv));
}

/**
 * Every non-excluded pair closer than the cutoff under the minimum image convention (found by scanning
 * lattice translations, independent of the code's sequential reduction) must be in the neighbor list.
 */
inline void checkNeighborCompleteness(int n, const std::vector<Vec3>& pos, const Vec3* box, bool periodic, double cutoff,
        const std::vector<std::set<int> >& exclusions, const std::vector<std::pair<unsigned int, unsigned int> >& list) {
    if (!enabled() || n > 400)
        return;
    std::set<std::pair<int, int> > listed;
    for (size_t k = 0; k < list.size(); k++) {
        int a = list[k].first, b = list[k].second;
        listed.insert(std::make_pair(a < b ? a : b, a < b ? b : a));
    }
    if (periodic) {
        // The minimum image is unique and found within one lattice shell only if each perpendicular cell width is at least twice the cutoff.
        Vec3 bc = box[1].cross(box[2]), ca = box[2].cross(box[0]), ab = box[0].cross(box[1]);
        double volume = std::fabs(box[0].dot(bc));
        if (volume/norm(bc) < 2*cutoff || volume/norm(ca) < 2*cutoff || volume/norm(ab) < 2*cutoff)
            return;
    }
    int range = periodic ? 2 : 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) {
            if (exclusions[i].count(j) || exclusions[j].count(i))
                continue;
            double best = 1e300;
            for (int a = -range; a <= range; a++)
                for (int b = -range; b <= range; b++)
                    for (int c = -range; c <= range; c++) {
                        Vec3 d = pos[j]-pos[i];
                        if (periodic)
                            d += box[0]*a + box[1]*b + box[2]*c;
                        double r = norm(d);
                        if (r < best)
                            best = r;
                    }
            if (best < cutoff*(1-1e-9) && !listed.count(std::make_pair(i, j))) {
                alarm("NBR-001", "neighbor_completeness", "missing pair at r " + num(best, cutoff));
                return;
            }
        }
}

/**
 * A particle with zero mass that is not a virtual site is not moved by an integration step.
 */
inline void checkMasslessFixed(const System& system, const std::vector<double>& masses, const std::vector<Vec3>& before, const std::vector<Vec3>& after) {
    // The masses are the ones the integrator uses (a later System::setParticleMass has no effect until reinitialization).
    if (!enabled() || before.size() != after.size() || (int) masses.size() != system.getNumParticles())
        return;
    for (int i = 0; i < system.getNumParticles(); i++)
        if (masses[i] == 0.0 && !system.isVirtualSite(i) && (after[i][0] != before[i][0] || after[i][1] != before[i][1] || after[i][2] != before[i][2])) {
            alarm("ZM-001", "massless_particle_fixed", "particle " + num(i) + " moved");
            return;
        }
}

/**
 * Total momentum of the massive particles vanishes after center-of-mass motion removal.
 */
inline void checkZeroMomentum(const std::vector<double>& masses, const std::vector<Vec3>& vel, double scaleBefore) {
    if (!enabled())
        return;
    Vec3 p;
    double scale = 0;
    for (size_t i = 0; i < masses.size(); i++)
        if (masses[i] != 0.0) {
            p += vel[i]*masses[i];
            scale += masses[i]*norm(vel[i]);
        }
    // The mean velocity is subtracted from velocities that may be much larger than the result, so rounding scales with the velocities before removal.
    scale += scaleBefore;
    if (norm(p) > 1024*EPS*(double) masses.size()*scale)
        alarm("CM-001", "momentum_removal", "residual momentum " + num(norm(p), scale));
}

inline Vec3 fractional(const Vec3* box, const Vec3& r) {
    double s2 = r[2]/box[2][2];
    double s1 = (r[1]-s2*box[2][1])/box[1][1];
    double s0 = (r[0]-s1*box[1][0]-s2*box[2][0])/box[0][0];
    return Vec3(s0, s1, s2);
}

/**
 * A barostat volume move is an affine deformation: the box has exactly the proposed volume, and the
 * geometric center of each independently scaled unit keeps its fractional coordinates in the box.
 */
inline void checkBarostatMove(double proposedVolume, const Vec3* oldBox, const Vec3* newBox, const std::vector<Vec3>& before,
        const std::vector<Vec3>& after, const std::vector<std::vector<int> >& units) {
    if (!enabled())
        return;
    double volume = newBox[0][0]*newBox[1][1]*newBox[2][2];
    if (std::fabs(volume-proposedVolume) > 64*EPS*proposedVolume)
        alarm("BARO-002", "barostat_volume", "box volume " + num(volume, proposedVolume));
    double worst = 0;
    // Fractional coordinates are formed by dividing absolute positions by the box (with shear terms), so their rounding
    // grows with how far the centers are from the origin in box units and with the tilt of the box.
    double mindiag = std::fmin(newBox[0][0], std::fmin(newBox[1][1], newBox[2][2]));
    double tilt = (std::fabs(newBox[1][0])+std::fabs(newBox[2][0])+std::fabs(newBox[2][1]))/mindiag;
    double reach = 0;
    for (size_t u = 0; u < units.size(); u++) {
        Vec3 c0, c1;
        for (size_t k = 0; k < units[u].size(); k++) {
            c0 += before[units[u][k]];
            c1 += after[units[u][k]];
        }
        c0 /= (double) units[u].size();
        c1 /= (double) units[u].size();
        Vec3 d = fractional(newBox, c1)-fractional(oldBox, c0);
        worst = std::fmax(worst, norm(d));
        reach = std::fmax(reach, std::fmax(norm(c0), norm(c1)));
    }
    if (worst > 1024*EPS*(double) before.size()*(1+reach/mindiag)*(1+tilt))
        alarm("BARO-003", "barostat_affine", "fractional shift " + num(worst));
}

/**
 * Running magnitude of the individual terms (pairs, bonds, torsions) a kernel added to the force array.  Rounding in the
 * accumulated sum scales with these, not with the (possibly cancelled) net force on each particle.
 */
struct TermScale {
    double force, torque;
};

inline TermScale& termScale() {
    static TermScale scale = {0.0, 0.0};
    return scale;
}

/**
 * Brackets one force term: on destruction adds the size of the change it made to the force on its particles.
 */
class TermProbe {
public:
    TermProbe(const std::vector<Vec3>& f, const std::vector<Vec3>& x, const int* idx, int n) : forces(f), pos(x), count(0) {
        if (!enabled())
            return;
        count = n < 8 ? n : 8;
        for (int k = 0; k < count; k++) {
            atoms[k] = idx[k];
            before[k] = f[idx[k]];
        }
    }
    ~TermProbe() {
        for (int k = 0; k < count; k++) {
            double d = norm(forces[atoms[k]]-before[k]);
            termScale().force += d;
            termScale().torque += norm(pos[atoms[k]])*d;
        }
    }
private:
    const std::vector<Vec3>& forces;
    const std::vector<Vec3>& pos;
    int count, atoms[8];
    Vec3 before[8];
};

/**
 * Records the force array before a kernel adds its contribution, then checks that the contribution is
 * translation invariant (sums to zero) and, when requested, rotation invariant (zero net torque).
 */
class ForceBalance {
public:
    ForceBalance(const std::vector<Vec3>& forces) {
        if (enabled()) {
            before = forces;
            termScale().force = termScale().torque = 0.0;
        }
    }
    void check(const char* forceId, const char* torqueId, const std::vector<Vec3>& after, const std::vector<Vec3>& pos, bool torque) const {
        if (!enabled() || before.size() != after.size())
            return;
        size_t n = after.size();
        Vec3 sum, tq;
        double fscale = 0, tscale = 0;
        for (size_t i = 0; i < n; i++) {
            Vec3 d = after[i]-before[i];
            sum += d;
            // The delta is formed by subtracting two accumulated arrays, so rounding scales with their magnitudes.
            double mag = norm(d) + norm(before[i]) + norm(after[i]);
            fscale += mag;
            if (torque) {
                tq += pos[i].cross(d);
                tscale += norm(pos[i])*mag;
            }
        }
        // Each pairwise/bonded term cancels exactly up to rounding; slack grows with the particle count and the force scale
        // (array magnitudes plus the magnitudes of the individual terms that were summed).
        fscale += termScale().force;
        tscale += termScale().torque;
        double ftol = 1024*EPS*(double) n*fscale;
        if (norm(sum) > ftol)
            alarm(forceId, "force_net_zero", "net force " + num(norm(sum), ftol));
        if (torque && norm(tq) > 1024*EPS*(double) n*tscale)
            alarm(torqueId, "force_net_torque_zero", "net torque " + num(norm(tq), 1024*EPS*(double) n*tscale));
    }
private:
    std::vector<Vec3> before;
};

} // namespace sanitizer
} // namespace OpenMM

#endif /*OPENMM_SANITIZER_H_*/
