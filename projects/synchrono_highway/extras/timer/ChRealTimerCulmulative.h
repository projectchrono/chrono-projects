

#ifndef CHREALTIMECULM_H
#define CHREALTIMECULM_H

#include <limits>

#include "chrono/core/ChTimer.h"

namespace chrono {

/// Class for a timer which attempts to enforce soft real-time.
class ChRealtimeCulmulative : public ChTimer<double> {
  public:
    /// Create the timer (outside the simulation loop, preferably just before beginning the loop)
    ChRealtimeCulmulative() { start(); }

    /// Call this function INSIDE the simulation loop, just ONCE per loop (preferably as the last call in the loop),
    /// passing it the integration step size used at this step. If the time elapsed over the last step (i.e., from
    /// the last call to Spin) is small than the integration step size, this function will spin in place until real time
    /// catches up with the simulation time, thus providing soft real-time capabilities.
    void Spin(double sim_time) {
        while (GetTimeSecondsIntermediate() < sim_time) {
        }
    }

    void Reset() {
        reset();
        start();
    }
};

}  // end namespace chrono

#endif
