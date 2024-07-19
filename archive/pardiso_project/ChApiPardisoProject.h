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
// Authors: Dario Mangoni, Radu Serban
// =============================================================================

#ifndef CHAPIPARDISOPROJECT_H
#define CHAPIPARDISOPROJECT_H

#include "chrono/ChVersion.h"
#include "chrono/core/ChPlatform.h"

// When compiling this library, remember to define CH_API_COMPILE_PARDISOPROJECT
// (so that the symbols with 'ChApiPardisoProject' in front of them will be
// marked as exported). Otherwise, just do not define it if you
// link the library to your code, and the symbols will be imported.

#if defined(CH_API_COMPILE_PARDISOPROJECT)
    #define ChApiPardisoProject ChApiEXPORT
#else
    #define ChApiPardisoProject ChApiIMPORT
#endif

/**
    @defgroup pardisoproject_module PardisoProject module
    @brief Module for the PardisoProject library direct solver

    Module provides access to the PardisoProject library. This library is
    currently used in Chrono for its parallel direct solver (Pardiso).

    For additional information, see:
    - the [installation guide](@ref module_pardisoproject_installation)
*/

#endif
