#ifndef DEMOS_UTILS_H
#define DEMOS_UTILS_H

#include <iostream>

#include "core/ChStream.h"
#include "chrono_multicore/physics/ChSystemMulticore.h"

// =============================================================================
// Utility function for displaying an ASCII progress bar for the quantity x
// which must be a value between 0 and n. The width 'w' represents the number
// of '=' characters corresponding to 100%.

static inline void progressbar(unsigned int x, unsigned int n, unsigned int w = 50) {
  if ((x != n) && (x % (n / 100 + 1) != 0))
    return;

  float ratio = x / (float)n;
  unsigned int c = (unsigned int)(ratio * w);

  std::cout << std::setw(3) << (int)(ratio * 100) << "% [";
  for (unsigned int x = 0; x < c; x++)
    std::cout << "=";
  for (unsigned int x = c; x < w; x++)
    std::cout << " ";
  std::cout << "]\r" << std::flush;
}

// =============================================================================
// Utility function to print to console a few important step statistics

static inline void TimingOutput(chrono::ChSystem* sys, chrono::ChStreamOutAsciiFile* ofile = NULL) {
  double TIME = sys->GetChTime();
  double STEP = sys->GetTimerStep();
  double BROD = sys->GetTimerCollisionBroad();
  double NARR = sys->GetTimerCollisionNarrow();
  double SOLVER = sys->GetTimerLSsolve();
  double UPDT = sys->GetTimerUpdate();
  double RESID = 0;
  int REQ_ITS = 0;
  int BODS = sys->GetNbodies();
  int CNTC = sys->GetNcontacts();
  if (chrono::ChSystemMulticore* mc_sys = dynamic_cast<chrono::ChSystemMulticore*>(sys)) {
      RESID = std::static_pointer_cast<chrono::ChIterativeSolverMulticore>(sys->GetSolver())->GetResidual();
      REQ_ITS = std::static_pointer_cast<chrono::ChIterativeSolverMulticore>(sys->GetSolver())->GetIterations();
      BODS = mc_sys->GetNbodies();
      CNTC = mc_sys->GetNcontacts();
  }

  if (ofile) {
      char buf[200];
      sprintf(buf, "%8.5f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7d  %7d  %7d  %7.4f\n", TIME, STEP, BROD, NARR, SOLVER,
              UPDT, BODS, CNTC, REQ_ITS, RESID);
      *ofile << buf;
  }

  printf("   %8.5f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7d | %7d | %7d | %7.4f\n", TIME, STEP, BROD, NARR, SOLVER,
         UPDT, BODS, CNTC, REQ_ITS, RESID);
}

#endif
