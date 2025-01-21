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

#include "chrono_multicore/solver/ChSolverMulticore.h"

using namespace chrono;

ChSchurProduct::ChSchurProduct() {
    data_manager = 0;
}
void ChSchurProduct::operator()(const DynamicVector<real>& x, DynamicVector<real>& output) {
    data_manager->system_timer.start("SchurProduct");

    const DynamicVector<real>& E = data_manager->host_data.E;

    uint num_rigid_contacts = data_manager->cd_data ? data_manager->cd_data->num_rigid_contacts : 0;
    uint num_unilaterals = data_manager->num_unilaterals;
    uint num_bilaterals = data_manager->num_bilaterals;
    output.reset();

    const CompressedMatrix<real>& D_T = data_manager->host_data.D_T;
    const CompressedMatrix<real>& Nschur = data_manager->host_data.Nschur;

    if (data_manager->settings.solver.local_solver_mode == data_manager->settings.solver.solver_mode) {
        if (data_manager->settings.solver.compute_N) {
            output = Nschur * x + E * x;
        } else {
            output = D_T * data_manager->host_data.M_invD * x + E * x;
        }

    } else {
        const SubMatrixType& D_n_T = _DNT_;
        const SubMatrixType& D_b_T = _DBT_;
        const SubMatrixType& M_invD_n = _MINVDN_;
        const SubMatrixType& M_invD_b = _MINVDB_;

        SubVectorType o_b = subvector(output, num_unilaterals, num_bilaterals);
        ConstSubVectorType x_b = subvector(x, num_unilaterals, num_bilaterals);
        ConstSubVectorType E_b = subvector(E, num_unilaterals, num_bilaterals);

        SubVectorType o_n = subvector(output, 0, num_rigid_contacts);
        ConstSubVectorType x_n = subvector(x, 0, num_rigid_contacts);
        ConstSubVectorType E_n = subvector(E, 0, num_rigid_contacts);

        switch (data_manager->settings.solver.local_solver_mode) {
            case SolverMode::BILATERAL: {
                o_b = D_b_T * (M_invD_b * x_b) + E_b * x_b;
            } break;

            case SolverMode::NORMAL: {
                blaze::DynamicVector<real> tmp = M_invD_b * x_b + M_invD_n * x_n;
                o_b = D_b_T * tmp + E_b * x_b;
                o_n = D_n_T * tmp + E_n * x_n;
            } break;

            case SolverMode::SLIDING: {
                const SubMatrixType& D_t_T = _DTT_;
                const SubMatrixType& M_invD_t = _MINVDT_;
                SubVectorType o_t = subvector(output, num_rigid_contacts, num_rigid_contacts * 2);
                ConstSubVectorType x_t = subvector(x, num_rigid_contacts, num_rigid_contacts * 2);
                ConstSubVectorType E_t = subvector(E, num_rigid_contacts, num_rigid_contacts * 2);

                blaze::DynamicVector<real> tmp = M_invD_b * x_b + M_invD_n * x_n + M_invD_t * x_t;
                o_b = D_b_T * tmp + E_b * x_b;
                o_n = D_n_T * tmp + E_n * x_n;
                o_t = D_t_T * tmp + E_t * x_t;

            } break;

            case SolverMode::SPINNING: {
                const SubMatrixType& D_t_T = _DTT_;
                const SubMatrixType& D_s_T = _DST_;
                const SubMatrixType& M_invD_t = _MINVDT_;
                const SubMatrixType& M_invD_s = _MINVDS_;
                SubVectorType o_t = subvector(output, num_rigid_contacts, num_rigid_contacts * 2);
                ConstSubVectorType x_t = subvector(x, num_rigid_contacts, num_rigid_contacts * 2);
                ConstSubVectorType E_t = subvector(E, num_rigid_contacts, num_rigid_contacts * 2);

                SubVectorType o_s = subvector(output, num_rigid_contacts * 3, num_rigid_contacts * 3);
                ConstSubVectorType x_s = subvector(x, num_rigid_contacts * 3, num_rigid_contacts * 3);
                ConstSubVectorType E_s = subvector(E, num_rigid_contacts * 3, num_rigid_contacts * 3);

                blaze::DynamicVector<real> tmp = M_invD_b * x_b + M_invD_n * x_n + M_invD_t * x_t + M_invD_s * x_s;
                o_b = D_b_T * tmp + E_b * x_b;
                o_n = D_n_T * tmp + E_n * x_n;
                o_t = D_t_T * tmp + E_t * x_t;
                o_s = D_s_T * tmp + E_s * x_s;

            } break;
        }
    }
    data_manager->system_timer.stop("SchurProduct");
}

void ChSchurProductBilateral::Setup(ChMulticoreDataManager* data_container_) {
    ChSchurProduct::Setup(data_container_);
    if (data_manager->num_bilaterals == 0) {
        return;
    }
    NschurB = _DBT_ * _MINVDB_;
}

void ChSchurProductBilateral::operator()(const DynamicVector<real>& x, DynamicVector<real>& output) {
    output = NschurB * x;
}
