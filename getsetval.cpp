#include "fem1dfunc.hpp"

void Neumann::setNuem(int x, int y, double z)
{
  elem = x;
  ndnum = y;
  val = z;
}

void Neumann::applyNuem(int **elemG, double *F)
{
  int jG;
  jG = elemG[elem][ndnum];
  F[jG] += val;
}

void Dirichlet::setDir(int x, double y)
{
  ndnum = x;
  val = y;  
}

void Dirichlet::applyDir(double **K, double *F, int n)
{

  for (int j = 0; j < n; j++)
        K[ndnum][j] = 0;

      K[ndnum][ndnum] = 1;
      F[ndnum] = 0;
      
      for (int i = 0; i < n; i++)
	{
	  if ( i != ndnum)
	    {
	      F[i] -= K[i][ndnum]*F[ndnum];
	      K[i][ndnum] = 0;
	    }
	}
      
}
