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
// Test for ChFunctionRecorder
// =============================================================================

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "chrono/utils/ChConstants.h"
#include "chrono/functions/ChFunctionInterp.h"

using namespace chrono;

double Reference(double x) {
    return std::cos(x * CH_C_2PI);
}

void Evaluate(const ChFunctionInterp& fun, double x) {
    double y = fun.GetVal(x);
    double yd = fun.GetDer(x);
    double ydd = fun.GetDer2(x);

    double y_ref = Reference(x);

    std::cout << "  " << x << "  " << y << "  " << yd << "  " << ydd;
    std::cout << "  " << y_ref << "  " << std::abs(y - y_ref) << std::endl;
}

int main(int argc, char* argv[]) {
    ChFunctionInterp fun;

    std::vector<double> x = {0, 0.1, 0.3, 0.32, 0.45, 0.5, 0.68, 0.7, 0.73, 0.79, 0.88, 1};

    std::cout << "Data\n";
    for (int i = 0; i < x.size(); i++) {
        std::cout << "  " << x[i] << "  " << Reference(x[i]) << std::endl;
    }

    // Insert data in recorder function in random order
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(x.begin(), x.end(), gen);
    for (int i = 0; i < x.size(); i++) {
        fun.AddPoint(x[i], Reference(x[i]));
    }

    // Check overwriting existing point
    double xx = 0.68 + 1e-14;
    fun.AddPoint(xx, Reference(xx), 1);

    ////Evaluate(fun, 0.03);
    ////Evaluate(fun, 0.02);

    // Evaluate below lower limit
    std::cout << "Below lower limit\n";
    Evaluate(fun, -0.1);

    // Evaluate above upper limit
    std::cout << "Above upper limit\n";
    Evaluate(fun, 1.1);

    // Evaluate in increasing order
    int n = 100;
    std::cout << "Increasing order\n";
    for (int i = 0; i <= n; i++) {
        double x = (i * 1.0) / n;
        Evaluate(fun, x);
    }

    ////Evaluate(fun, 0.85);
    ////Evaluate(fun, 0.2);
    ////Evaluate(fun, 0.08);

    // Evaluate in random order
    std::cout << "Random points\n";
    std::uniform_real_distribution<> dist(0, 1);
    for (int i = 0; i < n; i++) {
        Evaluate(fun, dist(gen));
    }
}
