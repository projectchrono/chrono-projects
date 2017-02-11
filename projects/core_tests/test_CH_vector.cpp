// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban
// =============================================================================
// Access test ChVector
// =============================================================================

#include <random>
#include <vector>

#include "chrono/core/ChVector.h"
#include "chrono/core/ChTimer.h"

using namespace chrono;

int main(int argc, char* argv[]) {
    ChTimer<> timer;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    size_t nv = 1 << 27;
    std::vector<ChVector<>> v;
    v.reserve(nv);

    for (size_t i = 0; i < nv; i++) {
        v.push_back(ChVector<>(dist(mt), dist(mt), dist(mt)));
        // v.push_back(ChVector<>(1.0, 2.0, -1.0));
    }

    std::cout << "number of vectors: " << v.size() << std::endl << std::endl;

    double xsum;
    double ysum;
    double zsum;

    xsum = 0;
    ysum = 0;
    zsum = 0;

    timer.reset();
    timer.start();
    for (auto& a : v) {
        xsum += a[0];
        ysum += a[1];
        zsum += a[2];
    }
    timer.stop();

    std::cout << "sums using operator[]: " << xsum << "  " << ysum << "  " << zsum << std::endl;
    std::cout << "time using operator[]: " << timer() << std::endl << std::endl;

    xsum = 0;
    ysum = 0;
    zsum = 0;

    timer.reset();
    timer.start();
    for (auto& a : v) {
        xsum += a.x();
        ysum += a.y();
        zsum += a.z();
    }
    timer.stop();

    std::cout << "sums using accessors: " << xsum << "  " << ysum << "  " << zsum << std::endl;
    std::cout << "time using accessors: " << timer() << std::endl << std::endl;

    ////xsum = 0;
    ////ysum = 0;
    ////zsum = 0;

    ////timer.reset();
    ////timer.start();
    ////for (auto& a : v) {
    ////    xsum += a.data[0];
    ////    ysum += a.data[1];
    ////    zsum += a.data[2];
    ////}
    ////timer.stop();

    ////std::cout << "sums using direct access: " << xsum << "  " << ysum << "  " << zsum << std::endl;
    ////std::cout << "time using direct access: " << timer() << std::endl << std::endl;
}
