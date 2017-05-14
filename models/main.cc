#include "TT06_RRG.h"

#include <fstream>


int main(int argc, char **argv)
{


  TT06_RRG t(1);

  double V = t.defaultVoltage();
  unsigned n = 10;
  double dt = 1.0 / n;

  std::ofstream file("out.txt");

  for(unsigned ii = 0; ii < 1000 * n; ii++) {

    double iStim = 0;

    if(ii * dt < 2) {
      iStim = 30;
    }

    double dV = t.computeRates(dt, V, 0, true);
    V += dV * dt + iStim * dt;

    if(ii % n == 0) {
      file << V << "\n";
    }

  }

  file.close();

}
