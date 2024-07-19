// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Dario Mangoni
// =============================================================================

#include "chrono_pardisoproject/ChSolverPardisoProject.h"
#include "chrono/utils/ChOpenMP.h"

namespace chrono {
ChSolverPardisoProject::ChSolverPardisoProject(int num_threads, ChPardisoProjectEngine::parproj_SYM symmetry)
    : m_engine(symmetry) {}

inline void ChSolverPardisoProject::SetMatrixSymmetryType(MatrixSymmetryType symmetry) {
    ChPardisoProjectEngine::parproj_SYM engine_symmetry;
    switch (symmetry) {
        case MatrixSymmetryType::GENERAL:
            engine_symmetry = ChPardisoProjectEngine::parproj_SYM::UNSYMMETRIC;
            break;
        case MatrixSymmetryType::STRUCTURAL_SYMMETRIC:
            engine_symmetry = ChPardisoProjectEngine::parproj_SYM::STRUCTURAL_SYMMETRIC;
            break;
        case MatrixSymmetryType::SYMMETRIC_INDEF:
            engine_symmetry = ChPardisoProjectEngine::parproj_SYM::SYMMETRIC_GENERAL;
            break;
        case MatrixSymmetryType::SYMMETRIC_POSDEF:
            engine_symmetry = ChPardisoProjectEngine::parproj_SYM::SYMMETRIC_POSDEF;
            break;
        default:
            std::cerr << "ChSolverPardisoProject does not support the matrix symmetry set." << std::endl
                      << "Rolling back to GENERAL" << std::endl;
            symmetry = MatrixSymmetryType::GENERAL;
            break;
    }
    m_symmetry = symmetry;
}

bool ChSolverPardisoProject::Setup(ChSystemDescriptor& sysd) {
    m_engine.SetZeroIndexedFormat();
    // TODO: yes, as ugly as it seems... the call to Setup is grabbed before the proper ChDirectSolverLS::Setup is
    // called so to have the chance to reset back the matrix to zero-based indexes we could avoid to set them back to
    // zero if we are sure that the sparsity pattern learner is disabled;
    return ChDirectSolverLS::Setup(sysd);
}

bool ChSolverPardisoProject::FactorizeMatrix() {
    m_engine.SetMatrix(m_mat);

    m_engine.PardisoProjectCall(ChPardisoProjectEngine::parproj_PHASE::ANALYZE_FACTORIZE);

    if (verbose) {
        if (m_engine.GetLastError() != 0) {
            std::cerr << std::endl << "ERROR during symbolic factorization: " << m_engine.GetLastError() << std::endl;
            exit(1);
        }
        std::endl << "Reordering completed ... " << std::endl;
        std::endl << "Number of nonzeros in factors  = " << m_engine.GetIPARM(17) << std::endl;
        std::endl << "Number of factorization MFLOPS = " << m_engine.GetIPARM(18) << std::endl;
    }

    return true;
}

bool ChSolverPardisoProject::SolveSystem() {
    m_engine.SetRhsVector(m_rhs);
    m_engine.SetSolutionVector(m_sol);
    m_engine.PardisoProjectCall(ChPardisoProjectEngine::parproj_PHASE::SOLVE);

    if (m_engine.GetLastError() != 0) {
        std::cout << "ERROR during solution: " << m_engine.GetLastError() << std::endl;
        exit(3);
    }
    return true;
}

void ChSolverPardisoProject::PrintErrorMessage() {
    std::cerr << "ERROR: " << m_engine.GetLastError() << std::endl;
}

}  // end of namespace chrono
