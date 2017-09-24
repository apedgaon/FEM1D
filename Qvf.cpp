#include "fem1dfunc.hpp"

void Qvf(int s, int *elem, double *node, double **B, double& Jac)
{
  double *N, *dNds;
  N = new double [2];
  dNds = new double[2];
  sfcn(s, N, dNds);

  Jac = 1.0/2*(node[elem[1]] - node[elem[0]]);
  double dNdx[2];
  for (int j = 0; j < 2; j++)
    dNdx[j] = dNds[j]*1/Jac;
  
  for (int j = 0; j < 2; j++)
    for (int k = 0; k < 2; k++)
      B[j][k] = (dNdx[j]*dNdx[k] - N[j]*N[k]);
  
  //Deallocate memory
  delete[] N;
  delete[] dNds;
}
