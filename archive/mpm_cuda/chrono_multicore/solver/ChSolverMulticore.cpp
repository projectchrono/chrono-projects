// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2016 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Hammad Mazhar
// =============================================================================
//
// This file contains the base class used for all multicore iterative solvers.
// All of the functions are defined here, with the implementation of each solver
// in it's specific cpp file.
//
// =============================================================================

#include "chrono_multicore/solver/ChSolverMulticore.h"

using namespace chrono;

void ChProjectConstraints::operator()(real* data) {
    data_manager->system_timer.start("ChSolverMulticore_Project");
    data_manager->rigid_rigid->Project(data);
    data_manager->node_container->Project(data);
    data_manager->system_timer.stop("ChSolverMulticore_Project");
}

ChSolverMulticore::ChSolverMulticore() {
    current_iteration = 0;
    rigid_rigid = NULL;
    three_dof = NULL;
    fem = NULL;
    bilateral = NULL;
}

//=================================================================================================================================

void ChSolverMulticore::ComputeSRhs(custom_vector<real>& gamma,
                                    const custom_vector<real>& rhs,
                                    custom_vector<real3>& vel_data,
                                    custom_vector<real3>& omg_data,
                                    custom_vector<real>& b) {
    // TODO change SHRS to use blaze
    // ComputeImpulses(gamma, vel_data, omg_data);
    // rigid_rigid->ComputeS(rhs, vel_data, omg_data, b);
}

bool init_eigen_vec = 0;

real ChSolverMulticore::LargestEigenValue(ChSchurProduct& SchurProduct, DynamicVector<real>& temp, real lambda) {
    eigen_vec.resize(temp.size());
    if (init_eigen_vec == 0) {
        eigen_vec = 1;
        init_eigen_vec = 1;
    }

    if (lambda != 0) {
        SchurProduct(eigen_vec, temp);
        eigen_vec = 1.0 / lambda * temp;
    }
    real lambda_old = 0;

    for (int i = 0; i < data_manager->settings.solver.max_power_iteration; i++) {
        SchurProduct(eigen_vec, temp);
        lambda = Sqrt((temp, temp));
        if (lambda == 0) {
            return 1;
        }
        printf("Lambda: %.20f \n", lambda);
        if (Abs(lambda_old - lambda) < data_manager->settings.solver.power_iter_tolerance) {
            break;
        }
        eigen_vec = 1.0 / lambda * temp;
        lambda_old = lambda;
    }
    return lambda;
}
