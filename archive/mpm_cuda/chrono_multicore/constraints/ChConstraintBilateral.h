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
// Authors: Hammad Mazhar, Radu Serban
// =============================================================================
//
// Handling of bilateral constraints for the system and Jacobian calculation
//
// =============================================================================

#pragma once

#include "chrono_multicore/ChDataManager.h"

namespace chrono {

/// @addtogroup multicore_constraint
/// @{

/// Bilateral (joint) constraints.
class CH_MULTICORE_API ChConstraintBilateral {
  public:
    ChConstraintBilateral() {}
    ~ChConstraintBilateral() {}

    void Setup(ChMulticoreDataManager* data_container_) { data_manager = data_container_; }

    /// Compute the vector of corrections.
    void Build_b();
    /// Compute the diagonal compliance matrix.
    void Build_E();
    /// Compute the jacobian matrix, no allocation is performed here,
    /// GenerateSparsity should take care of that.
    void Build_D();

    //// Fill-in the non zero entries in the bilateral jacobian with ones.
    // This operation is sequential.
    void GenerateSparsity();

    ChMulticoreDataManager* data_manager;  ///< Pointer to the system's data manager.
};

/// @} multicore_constraint

}  // end namespace chrono
