/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2014-2019 Stanford University and the Authors.      *
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
#ifdef _MSC_VER
    // Prevent Windows from defining macros that interfere with other code.
    #define NOMINMAX
#endif

#include <algorithm>

#include "ReferenceTabulatedFunction.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/SplineFitter.h"

#include <cmath>
#include "openmm/internal/Sanitizer.h"

#ifdef _MSC_VER
#if _MSC_VER < 1800
/**
 * We need to define this ourselves, since Visual Studio is missing round() from cmath.
 */
static int round(double x) {
    return (int) (x+0.5);
}
#endif
#endif

static double wrap(double t, double min, double max) {
    double L = max - min;
    double s = (t - min)/L;
    return min + L*(s - floor(s));
}

using namespace OpenMM;
using namespace std;
using Lepton::CustomFunction;


// A tabulated function interpolates its data: evaluating it at a grid node returns the stored value.  For a periodic
// dimension the last node is identified with the first and is not checked.
static void checkNodes(const char* what, const Lepton::CustomFunction& f, const std::vector<double>& values, const std::vector<std::vector<double> >& limits, bool periodic) {
    // limits[d] = {min, max, count}; node positions are recomputed here rather than taken from the implementation.
    std::vector<std::vector<double> > nodes(limits.size());
    for (size_t d = 0; d < limits.size(); d++)
        for (int i = 0; i < (int) limits[d][2]; i++)
            nodes[d].push_back(limits[d][0] + i*(limits[d][1]-limits[d][0])/(limits[d][2]-1));
    if (!OpenMM::sanitizer::enabled() || values.size() > 100000)
        return;
    double scale = 0;
    for (double v : values)
        scale = std::max(scale, fabs(v));
    int dims = nodes.size();
    std::vector<int> n(dims), idx(dims, 0);
    size_t total = 1;
    for (int d = 0; d < dims; d++) {
        n[d] = nodes[d].size() - (periodic ? 1 : 0);
        total *= n[d];
    }
    for (size_t count = 0; count < total; count++) {
        double args[3];
        size_t flat = 0, stride = 1;
        for (int d = 0; d < dims; d++) {
            args[d] = nodes[d][idx[d]];
            flat += idx[d]*stride;
            stride *= nodes[d].size();
        }
        double got = f.evaluate(args), want = values[flat];
        if (fabs(got-want) > 1e3*OpenMM::sanitizer::EPS*(scale+1e-300)) {
            OpenMM::sanitizer::alarm("TAB-001", "interpolation_at_nodes", std::string(what) + " node value " + OpenMM::sanitizer::num(got, want));
            return;
        }
        for (int d = 0; d < dims; d++) {
            if (++idx[d] < n[d])
                break;
            idx[d] = 0;
        }
    }
}

extern "C" OPENMM_EXPORT CustomFunction* createReferenceTabulatedFunction(const TabulatedFunction& function) {
    CustomFunction* fn;
    if (dynamic_cast<const Continuous1DFunction*>(&function) != NULL)
        fn = new ReferenceContinuous1DFunction(dynamic_cast<const Continuous1DFunction&>(function));
    else if (dynamic_cast<const Continuous2DFunction*>(&function) != NULL)
        fn = new ReferenceContinuous2DFunction(dynamic_cast<const Continuous2DFunction&>(function));
    else if (dynamic_cast<const Continuous3DFunction*>(&function) != NULL)
        fn = new ReferenceContinuous3DFunction(dynamic_cast<const Continuous3DFunction&>(function));
    else if (dynamic_cast<const Discrete1DFunction*>(&function) != NULL)
        fn = new ReferenceDiscrete1DFunction(dynamic_cast<const Discrete1DFunction&>(function));
    else if (dynamic_cast<const Discrete2DFunction*>(&function) != NULL)
        fn = new ReferenceDiscrete2DFunction(dynamic_cast<const Discrete2DFunction&>(function));
    else if (dynamic_cast<const Discrete3DFunction*>(&function) != NULL)
        fn = new ReferenceDiscrete3DFunction(dynamic_cast<const Discrete3DFunction&>(function));
    else
        throw OpenMMException("createReferenceTabulatedFunction: Unknown function type");
    return new SharedFunctionWrapper(shared_ptr<const CustomFunction>(fn));
}

ReferenceContinuous1DFunction::ReferenceContinuous1DFunction(const Continuous1DFunction& function) : function(function) {
    periodic = function.getPeriodic();
    function.getFunctionParameters(values, min, max);
    int numValues = values.size();
    x.resize(numValues);
    for (int i = 0; i < numValues; i++)
        x[i] = min+i*(max-min)/(numValues-1);
    SplineFitter::createSpline(x, values, periodic, derivs);
    {
        std::vector<std::vector<double> > lim; lim.push_back({min, max, (double) values.size()});
        checkNodes("1D", *this, values, lim, periodic);
    }
}

ReferenceContinuous1DFunction::ReferenceContinuous1DFunction(const ReferenceContinuous1DFunction& other) : function(other.function) {
    periodic = function.getPeriodic();
    function.getFunctionParameters(values, min, max);
    x = other.x;
    values = other.values;
    derivs = other.derivs;
}

int ReferenceContinuous1DFunction::getNumArguments() const {
    return 1;
}

double ReferenceContinuous1DFunction::evaluate(const double* arguments) const {
    double t = periodic ? wrap(arguments[0], min, max) : arguments[0];
    if (t < min || t > max)
        return 0.0;
    return SplineFitter::evaluateSpline(x, values, derivs, t);
}

double ReferenceContinuous1DFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    double t = periodic ? wrap(arguments[0], min, max) : arguments[0];
    if (t < min || t > max)
        return 0.0;
    return SplineFitter::evaluateSplineDerivative(x, values, derivs, t);
}

CustomFunction* ReferenceContinuous1DFunction::clone() const {
    return new ReferenceContinuous1DFunction(*this);
}

ReferenceContinuous2DFunction::ReferenceContinuous2DFunction(const Continuous2DFunction& function) : function(function) {
    periodic = function.getPeriodic();
    function.getFunctionParameters(xsize, ysize, values, xmin, xmax, ymin, ymax);
    x.resize(xsize);
    y.resize(ysize);
    for (int i = 0; i < xsize; i++)
        x[i] = xmin+i*(xmax-xmin)/(xsize-1);
    for (int i = 0; i < ysize; i++)
        y[i] = ymin+i*(ymax-ymin)/(ysize-1);
    SplineFitter::create2DSpline(x, y, values, periodic, c);
    {
        std::vector<std::vector<double> > lim; lim.push_back({xmin, xmax, (double) xsize}); lim.push_back({ymin, ymax, (double) ysize});
        checkNodes("2D", *this, values, lim, periodic);
    }
}

ReferenceContinuous2DFunction::ReferenceContinuous2DFunction(const ReferenceContinuous2DFunction& other) : function(other.function) {
    periodic = function.getPeriodic();
    function.getFunctionParameters(xsize, ysize, values, xmin, xmax, ymin, ymax);
    x = other.x;
    y = other.y;
    values = other.values;
    c = other.c;
}

int ReferenceContinuous2DFunction::getNumArguments() const {
    return 2;
}

double ReferenceContinuous2DFunction::evaluate(const double* arguments) const {
    double u = periodic ? wrap(arguments[0], xmin, xmax) : arguments[0];
    if (u < xmin || u > xmax)
        return 0.0;
    double v = periodic ? wrap(arguments[1], ymin, ymax) : arguments[1];
    if (v < ymin || v > ymax)
        return 0.0;
    return SplineFitter::evaluate2DSpline(x, y, values, c, u, v);
}

double ReferenceContinuous2DFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    double u = periodic ? wrap(arguments[0], xmin, xmax) : arguments[0];
    if (u < xmin || u > xmax)
        return 0.0;
    double v = periodic ? wrap(arguments[1], ymin, ymax) : arguments[1];
    if (v < ymin || v > ymax)
        return 0.0;
    double dx, dy;
    SplineFitter::evaluate2DSplineDerivatives(x, y, values, c, u, v, dx, dy);
    if (derivOrder[0] == 1 && derivOrder[1] == 0)
        return dx;
    if (derivOrder[0] == 0 && derivOrder[1] == 1)
        return dy;
    throw OpenMMException("ReferenceContinuous2DFunction: Unsupported derivative order");
}

CustomFunction* ReferenceContinuous2DFunction::clone() const {
    return new ReferenceContinuous2DFunction(*this);
}

ReferenceContinuous3DFunction::ReferenceContinuous3DFunction(const Continuous3DFunction& function) : function(function) {
    periodic = function.getPeriodic();
    function.getFunctionParameters(xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax);
    x.resize(xsize);
    y.resize(ysize);
    z.resize(zsize);
    for (int i = 0; i < xsize; i++)
        x[i] = xmin+i*(xmax-xmin)/(xsize-1);
    for (int i = 0; i < ysize; i++)
        y[i] = ymin+i*(ymax-ymin)/(ysize-1);
    for (int i = 0; i < zsize; i++)
        z[i] = zmin+i*(zmax-zmin)/(zsize-1);
    SplineFitter::create3DSpline(x, y, z, values, periodic, c);
    {
        std::vector<std::vector<double> > lim; lim.push_back({xmin, xmax, (double) xsize}); lim.push_back({ymin, ymax, (double) ysize}); lim.push_back({zmin, zmax, (double) zsize});
        checkNodes("3D", *this, values, lim, periodic);
    }
}

ReferenceContinuous3DFunction::ReferenceContinuous3DFunction(const ReferenceContinuous3DFunction& other) : function(other.function) {
    periodic = function.getPeriodic();
    function.getFunctionParameters(xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax);
    x = other.x;
    y = other.y;
    z = other.z;
    values = other.values;
    c = other.c;
}

int ReferenceContinuous3DFunction::getNumArguments() const {
    return 3;
}

double ReferenceContinuous3DFunction::evaluate(const double* arguments) const {
    double u = periodic ? wrap(arguments[0], xmin, xmax) : arguments[0];
    if (u < xmin || u > xmax)
        return 0.0;
    double v = periodic ? wrap(arguments[1], ymin, ymax) : arguments[1];
    if (v < ymin || v > ymax)
        return 0.0;
    double w = periodic ? wrap(arguments[2], zmin, zmax) : arguments[2];
    if (w < zmin || w > zmax)
        return 0.0;
    return SplineFitter::evaluate3DSpline(x, y, z, values, c, u, v, w);
}

double ReferenceContinuous3DFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    double u = periodic ? wrap(arguments[0], xmin, xmax) : arguments[0];
    if (u < xmin || u > xmax)
        return 0.0;
    double v = periodic ? wrap(arguments[1], ymin, ymax) : arguments[1];
    if (v < ymin || v > ymax)
        return 0.0;
    double w = periodic ? wrap(arguments[2], zmin, zmax) : arguments[2];
    if (w < zmin || w > zmax)
        return 0.0;
    double dx, dy, dz;
    SplineFitter::evaluate3DSplineDerivatives(x, y, z, values, c, u, v, w, dx, dy, dz);
    if (derivOrder[0] == 1 && derivOrder[1] == 0 && derivOrder[2] == 0)
        return dx;
    if (derivOrder[0] == 0 && derivOrder[1] == 1 && derivOrder[2] == 0)
        return dy;
    if (derivOrder[0] == 0 && derivOrder[1] == 0 && derivOrder[2] == 1)
        return dz;
    throw OpenMMException("ReferenceContinuous3DFunction: Unsupported derivative order");
}

CustomFunction* ReferenceContinuous3DFunction::clone() const {
    return new ReferenceContinuous3DFunction(*this);
}

ReferenceDiscrete1DFunction::ReferenceDiscrete1DFunction(const Discrete1DFunction& function) : function(function) {
    function.getFunctionParameters(values);
}

int ReferenceDiscrete1DFunction::getNumArguments() const {
    return 1;
}

double ReferenceDiscrete1DFunction::evaluate(const double* arguments) const {
    int i = (int) round(arguments[0]);
    i = max(0, min((int) values.size()-1, i));
    return values[i];
}

double ReferenceDiscrete1DFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    return 0.0;
}

CustomFunction* ReferenceDiscrete1DFunction::clone() const {
    return new ReferenceDiscrete1DFunction(function);
}

ReferenceDiscrete2DFunction::ReferenceDiscrete2DFunction(const Discrete2DFunction& function) : function(function) {
    function.getFunctionParameters(xsize, ysize, values);
}

int ReferenceDiscrete2DFunction::getNumArguments() const {
    return 2;
}

double ReferenceDiscrete2DFunction::evaluate(const double* arguments) const {
    int i = (int) round(arguments[0]);
    int j = (int) round(arguments[1]);
    i = max(0, min(xsize-1, i));
    j = max(0, min(ysize-1, j));
    return values[i+j*xsize];
}

double ReferenceDiscrete2DFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    return 0.0;
}

CustomFunction* ReferenceDiscrete2DFunction::clone() const {
    return new ReferenceDiscrete2DFunction(function);
}

ReferenceDiscrete3DFunction::ReferenceDiscrete3DFunction(const Discrete3DFunction& function) : function(function) {
    function.getFunctionParameters(xsize, ysize, zsize, values);
}

int ReferenceDiscrete3DFunction::getNumArguments() const {
    return 3;
}

double ReferenceDiscrete3DFunction::evaluate(const double* arguments) const {
    int i = (int) round(arguments[0]);
    int j = (int) round(arguments[1]);
    int k = (int) round(arguments[2]);
    i = max(0, min(xsize-1, i));
    j = max(0, min(ysize-1, j));
    k = max(0, min(zsize-1, k));
    return values[i+(j+k*ysize)*xsize];
}

double ReferenceDiscrete3DFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    return 0.0;
}

CustomFunction* ReferenceDiscrete3DFunction::clone() const {
    return new ReferenceDiscrete3DFunction(function);
}

SharedFunctionWrapper::SharedFunctionWrapper(shared_ptr<const CustomFunction> pointer) : pointer(pointer) {
}

int SharedFunctionWrapper::getNumArguments() const {
    return pointer->getNumArguments();
}

double SharedFunctionWrapper::evaluate(const double* arguments) const {
    return pointer->evaluate(arguments);
}

double SharedFunctionWrapper::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    return pointer->evaluateDerivative(arguments, derivOrder);
}

CustomFunction* SharedFunctionWrapper::clone() const {
    return new SharedFunctionWrapper(pointer);
}
