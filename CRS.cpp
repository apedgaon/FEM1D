#include <cmath>
#include "fem1dfunc.hpp"

void CRS (double **A, double *values, int *col_in, int *row_p, int n, int *cs)
{
  for (int i=0; i<n; i++)
    {
      row_p[i] = *cs;
      for (int j=0; j<n; j++)
	{
	  if (fabs(A[i][j] - 0.0) > 1e-5)
	    {
	      values[*cs] = A[i][j];
	      col_in[*cs] = j;
	      *cs = *cs + 1;
	    }
	}      
    }
  row_p[n] = *cs;
}
