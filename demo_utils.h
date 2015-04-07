#ifndef DEMOS_UTILS_H
#define DEMOS_UTILS_H

#include <iostream>


// =============================================================================
// Utility function for displaying an ASCII progress bar for the quantity x
// which must be a value between 0 and n. The width 'w' represents the number
// of '=' characters corresponding to 100%.

static inline void progressbar(unsigned int x, unsigned int n, unsigned int w = 50) {
  if ((x != n) && (x % (n / 100 + 1) != 0))
    return;

  float ratio = x / (float)n;
  int c = ratio * w;

  std::cout << std::setw(3) << (int)(ratio * 100) << "% [";
  for (int x = 0; x < c; x++)
    std::cout << "=";
  for (int x = c; x < w; x++)
    std::cout << " ";
  std::cout << "]\r" << std::flush;
}



#endif
