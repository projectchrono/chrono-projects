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
//
// Performance test for MatrMultiplyAVX and MatrMultiplyTAVX.
//
// =============================================================================

#include "chrono/core/ChMatrixDynamic.h"
#include "chrono/core/ChTimer.h"

#include "chrono_thirdparty/filesystem/path.h"

#include "../BaseTest.h"

using namespace chrono;

// -----------------------------------------------------------------------------

// Test class
class MatMultAVX : public BaseTest {
  public:
    MatMultAVX(const std::string& testName,         // name of current test
               const std::string& testProjectName,  // name of tested module
               bool transposed                      // if true, test A * B^T
               )
        : BaseTest(testName, testProjectName), m_execTime(0), m_transposed(transposed) {}

    ~MatMultAVX() {}

    // Override corresponding functions in BaseTest
    virtual bool execute() override;
    virtual double getExecutionTime() const override { return m_execTime; }

private:

    bool CheckMatMult(int M, int N, int K);
    bool CheckMatMultT(int M, int N, int K);

    double m_execTime;
    bool m_transposed;

    static const int m_num_runs;
    static const double m_tolerance;
};

const int MatMultAVX::m_num_runs = 10000;
const double MatMultAVX::m_tolerance = 1e-12;

bool MatMultAVX::execute() {
    bool passed = true;
    ChTimer<double> timer;

    if (m_transposed) {
        timer.start();
        for (int run = 0; run < m_num_runs; run++) {
            for (int N = 20; N <= 23; N++) {
                for (int M = 8; M <= 11; M++) {
                    passed &= CheckMatMultT(M, N, 24);
                }
            }
        }
        timer.stop();
        m_execTime = timer.GetTimeSeconds();
    } else {
        timer.start();
        for (int run = 0; run < m_num_runs; run++) {
            for (int N = 20; N <= 23; N++) {
                for (int M = 8; M <= 11; M++) {
                    passed &= CheckMatMult(M, N, 24);
                }
            }
        }
        timer.stop();
        m_execTime = timer.GetTimeSeconds();
    }

    return passed;
}

// -----------------------------------------------------------------------------

// Check multiplication A*B of random matrices A (MxN) and B (NxK)
bool MatMultAVX::CheckMatMult(int M, int N, int K) {
    ChMatrixDynamic<double> A(M, N);
    ChMatrixDynamic<double> B(N, K);

    A.FillRandom(10, -10);
    B.FillRandom(10, -10);

    ChMatrixDynamic<double> ref(M, K);
    ref.MatrMultiply(A, B);

    ChMatrixDynamic<double> avx(M, K);
    avx.MatrMultiplyAVX(A, B);

    return avx.Equals(ref, m_tolerance);
}

// Check multiplication A*B' of random matrices A (MxN) and B (KxN)
bool MatMultAVX::CheckMatMultT(int M, int N, int K) {
    ChMatrixDynamic<double> A(M, N);
    ChMatrixDynamic<double> B(K, N);

    A.FillRandom(10, -10);
    B.FillRandom(10, -10);

    ChMatrixDynamic<double> ref(M, K);
    ref.MatrMultiplyT(A, B);

    ChMatrixDynamic<double> avx(M, K);
    avx.MatrMultiplyTAVX(A, B);

    return avx.Equals(ref, m_tolerance);
}

int main(int argc, char* argv[]) {
    std::string out_dir = "../METRICS";
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return 1;
    }

    MatMultAVX mult("metrics_CH_matMultAVX", "Chrono::Engine", false);
    MatMultAVX multT("metrics_CH_matMultAVX_T", "Chrono::Engine", true);

    mult.setOutDir(out_dir);
    mult.setVerbose(true);
    mult.run();
    mult.print();

    multT.setOutDir(out_dir);
    multT.setVerbose(true);
    multT.run();
    multT.print();

    return 0;
}
