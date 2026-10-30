/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2013-2015 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "ReferenceConstraints.h"
#include "ReferenceCCMAAlgorithm.h"
#include "ReferenceSETTLEAlgorithm.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/OpenMMException.h"
#include <map>
#include <utility>
#include "openmm/internal/Sanitizer.h"
#include <vector>

using namespace OpenMM;
using namespace std;

ReferenceConstraints::ReferenceConstraints(const System& system) : ccma(NULL), settle(NULL) {
    int numParticles = system.getNumParticles();
    vector<double> masses(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);

    // Record the set of constraints and how many constraints each atom is involved in.

    vector<int> atom1;
    vector<int> atom2;
    vector<double> distance;
    vector<int> constraintCount(numParticles, 0);
    for (int i = 0; i < system.getNumConstraints(); i++) {
        int p1, p2;
        double d;
        system.getConstraintParameters(i, p1, p2, d);
        if (masses[p1] != 0 || masses[p2] != 0) {
            atom1.push_back(p1);
            atom2.push_back(p2);
            distance.push_back(d);
            constraintCount[p1]++;
            constraintCount[p2]++;
        }
    }

    sanAtom1 = atom1;
    sanAtom2 = atom2;
    sanDistance = distance;

    // Identify clusters of three atoms that can be treated with SETTLE.  First, for every
    // atom that might be part of such a cluster, make a list of the two other atoms it is
    // connected to.

    vector<map<int, float> > settleConstraints(numParticles);
    for (int i = 0; i < (int)atom1.size(); i++) {
        if (constraintCount[atom1[i]] == 2 && constraintCount[atom2[i]] == 2) {
            settleConstraints[atom1[i]][atom2[i]] = (float) distance[i];
            settleConstraints[atom2[i]][atom1[i]] = (float) distance[i];
        }
    }

    // Now remove the ones that don't actually form closed loops of three atoms.

    vector<int> settleClusters;
    for (int i = 0; i < (int)settleConstraints.size(); i++) {
        if (settleConstraints[i].size() == 2) {
            int partner1 = settleConstraints[i].begin()->first;
            int partner2 = (++settleConstraints[i].begin())->first;
            if (settleConstraints[partner1].size() != 2 || settleConstraints[partner2].size() != 2 ||
                    settleConstraints[partner1].find(partner2) == settleConstraints[partner1].end())
                settleConstraints[i].clear();
            else if (i < partner1 && i < partner2)
                settleClusters.push_back(i);
        }
        else
            settleConstraints[i].clear();
    }

    // Record the SETTLE clusters.

    vector<bool> isSettleAtom(numParticles, false);
    if (settleClusters.size() > 0) {
        vector<int> atom1;
        vector<int> atom2;
        vector<int> atom3;
        vector<double> distance1;
        vector<double> distance2;
        for (int i = 0; i < settleClusters.size(); i++) {
            int p1 = settleClusters[i];
            int p2 = settleConstraints[p1].begin()->first;
            int p3 = (++settleConstraints[p1].begin())->first;
            float dist12 = settleConstraints[p1].find(p2)->second;
            float dist13 = settleConstraints[p1].find(p3)->second;
            float dist23 = settleConstraints[p2].find(p3)->second;
            if (dist12 == dist13) {
                // p1 is the central atom
                atom1.push_back(p1);
                atom2.push_back(p2);
                atom3.push_back(p3);
                distance1.push_back(dist12);
                distance2.push_back(dist23);
            }
            else if (dist12 == dist23) {
                // p2 is the central atom
                atom1.push_back(p2);
                atom2.push_back(p1);
                atom3.push_back(p3);
                distance1.push_back(dist12);
                distance2.push_back(dist13);
            }
            else if (dist13 == dist23) {
                // p3 is the central atom
                atom1.push_back(p3);
                atom2.push_back(p1);
                atom3.push_back(p2);
                distance1.push_back(dist13);
                distance2.push_back(dist12);
            }
            else
                continue; // We can't handle this with SETTLE
            isSettleAtom[p1] = true;
            isSettleAtom[p2] = true;
            isSettleAtom[p3] = true;
        }
        if (atom1.size() > 0)
            settle = new ReferenceSETTLEAlgorithm(atom1, atom2, atom3, distance1, distance2, masses);
    }

    // All other constraints are handled with CCMA.

    vector<int> ccmaConstraints;
    for (unsigned i = 0; i < atom1.size(); i++)
        if (!isSettleAtom[atom1[i]])
            ccmaConstraints.push_back(i);
    int numCCMA = ccmaConstraints.size();
    if (numCCMA > 0) {
        // Record particles and distances for CCMA.
        
        vector<pair<int, int> > ccmaIndices(numCCMA);
        vector<double> ccmaDistance(numCCMA);
        for (int i = 0; i < numCCMA; i++) {
            int index = ccmaConstraints[i];
            ccmaIndices[i] = make_pair(atom1[index], atom2[index]);
            ccmaDistance[i] = distance[index];
        }

        // Look up angles for CCMA.
        
        vector<ReferenceCCMAAlgorithm::AngleInfo> angles;
        for (int i = 0; i < system.getNumForces(); i++) {
            const HarmonicAngleForce* force = dynamic_cast<const HarmonicAngleForce*>(&system.getForce(i));
            if (force != NULL) {
                for (int j = 0; j < force->getNumAngles(); j++) {
                    int atom1, atom2, atom3;
                    double angle, k;
                    force->getAngleParameters(j, atom1, atom2, atom3, angle, k);
                    angles.push_back(ReferenceCCMAAlgorithm::AngleInfo(atom1, atom2, atom3, angle));
                }
            }
        }
        
        // Create the CCMA object.
        
        ccma = new ReferenceCCMAAlgorithm(numParticles, numCCMA, ccmaIndices, ccmaDistance, masses, angles, 0.02);
    }
}

ReferenceConstraints::~ReferenceConstraints() {
    if (ccma != NULL)
        delete ccma;
    if (settle != NULL)
        delete settle;
}

void ReferenceConstraints::apply(vector<OpenMM::Vec3>& atomCoordinates, vector<OpenMM::Vec3>& atomCoordinatesP, vector<double>& inverseMasses, double tolerance) {
    vector<Vec3> before;
    if (sanitizer::enabled())
        before = atomCoordinatesP;
    if (ccma != NULL)
        ccma->apply(atomCoordinates, atomCoordinatesP, inverseMasses, tolerance);
    if (settle != NULL)
        settle->apply(atomCoordinates, atomCoordinatesP, inverseMasses, tolerance);
    if (!sanitizer::enabled())
        return;

    // The solver stops when every constrained distance is within the tolerance.
    double worst = 0;
    for (size_t i = 0; i < sanAtom1.size(); i++) {
        Vec3 d = atomCoordinatesP[sanAtom1[i]]-atomCoordinatesP[sanAtom2[i]];
        double err = fabs(sanitizer::norm(d)-sanDistance[i])/sanDistance[i];
        worst = max(worst, err);
    }
    if (worst > 2*tolerance + 1e-12)
        sanitizer::alarm("CON-001", "constraint_postcondition", "relative distance error " + sanitizer::num(worst, 2*tolerance));

    // Constraint forces are internal: the mass-weighted displacement they cause sums to zero.
    Vec3 sum;
    double scale = 0;
    for (size_t i = 0; i < atomCoordinatesP.size(); i++)
        if (inverseMasses[i] != 0.0) {
            Vec3 dx = atomCoordinatesP[i]-before[i];
            sum += dx/inverseMasses[i];
            // Positions are rebuilt from absolute coordinates, so rounding scales with |x|, not only with the displacement.
            scale += (sanitizer::norm(dx)+sanitizer::norm(atomCoordinatesP[i]))/inverseMasses[i];
        }
    if (sanitizer::norm(sum) > 1024*sanitizer::EPS*(double) atomCoordinatesP.size()*scale)
        sanitizer::alarm("CON-003", "constraint_momentum", "net displacement impulse " + sanitizer::num(sanitizer::norm(sum), scale));
}

void ReferenceConstraints::applyToVelocities(vector<OpenMM::Vec3>& atomCoordinates, vector<OpenMM::Vec3>& velocities, vector<double>& inverseMasses, double tolerance) {
    if (ccma != NULL)
        ccma->applyToVelocities(atomCoordinates, velocities, inverseMasses, tolerance);
    if (settle != NULL)
        settle->applyToVelocities(atomCoordinates, velocities, inverseMasses, tolerance);
    if (!sanitizer::enabled())
        return;

    // After projection no constrained pair moves along its bond.
    double vmax = 0;
    for (size_t i = 0; i < velocities.size(); i++)
        vmax = max(vmax, sanitizer::norm(velocities[i]));
    double worst = 0;
    for (size_t i = 0; i < sanAtom1.size(); i++) {
        Vec3 e = atomCoordinates[sanAtom1[i]]-atomCoordinates[sanAtom2[i]];
        double len = sanitizer::norm(e);
        if (len == 0)
            continue;
        Vec3 dv = velocities[sanAtom1[i]]-velocities[sanAtom2[i]];
        worst = max(worst, fabs(dv.dot(e))/len);
    }
    if (worst > 2*tolerance*vmax + 1e3*sanitizer::EPS*vmax)
        sanitizer::alarm("CON-002", "constraint_postcondition", "relative speed along bond " + sanitizer::num(worst, vmax));
}
