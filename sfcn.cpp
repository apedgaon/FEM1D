#include "fem1dfunc.hpp"

void sfcn(double s, double *N, double *dNds)
{
  N[0] = (1 - s)/2;
  N[1] = (1 + s)/2;

  dNds[0] = -0.5;
  dNds[1] = 0.5;
}
