#include <iostream>
#include <cmath>
#include "fem1dfunc.hpp"

void CG_Solve(double **A, double *b, double*u, int n)
{

  double *res, *p, *zres;
  res = new double [n];
  p = new double [n];
  zres = new double [n];
  double norm_res = 0.0;
  double norm_b = 0.0;
  
  double **M_inv;
  M_inv = new double *[n];
  
  for (int i=0; i<n; i++)
    {
      M_inv[i] = new double [n];
    }
  
  for (int i=0; i<n; i++)
    {
      for (int j=0; j<n; j++)
	{
	  if (i == j)
	    {
	      M_inv[i][j] = 1.0/A[i][j];
	    }
	  else
	    {
	      M_inv[i][j] = 0.0;
	    }
	}
    }

  double *values;
  int *col_in, *row_p;
  int cs = 0;
  values = new double [n*n];
  col_in = new int [n*n];
  row_p = new int [n+1];  
  
  CRS(A, values, col_in, row_p, n, &cs);  

  for (int i=0; i<n; i++)
    {
      res[i] = b[i];
      zres[i] = M_inv[i][i]*res[i];
      p[i] = zres[i];
      norm_res += pow(res[i],2);      
    }
  
  norm_res = sqrt(norm_res);
  norm_b = norm_res;
  double eps = 1e-6;
  int k = 0;
  double alpha, num, numn, *denv, den, Beta;
  denv = new double [n];

  while (k < 1e5 && norm_res > eps*norm_b)
    {
      num = 0.0;
      den = 0.0;
      
      for (int i=0; i<n; i++)
	{
	  denv[i] = 0.0;
	}

      for (int i=0; i<n; i++)
	{
	  num += res[i]*zres[i];
	  
	  for (int j=row_p[i]; j<row_p[i+1]; j++)
	    {
	      denv[i] += values[j]*p[col_in[j]];
	    }
	  den += p[i]*denv[i];
	}
      
      alpha = num/den;
      numn = 0.0;

      for (int i=0; i<n; i++)
	{
	  u[i] = u[i] + alpha*p[i];
	  res[i] = res[i] - alpha*denv[i];
	  zres[i] = M_inv[i][i]*res[i];
	  numn += res[i]*zres[i];
	  
	}
      
      Beta = numn/num;
      norm_res = 0.0;

      for (int i=0; i<n; i++)
	{
	  p[i] = zres[i] + Beta*p[i];
	  norm_res += res[i]*res[i];
	}
      
      norm_res = sqrt(norm_res);
      k++;      
    }

  std::cout << "iterations = " << k << "\n";

  for (int i=0; i<n; i++)
    {
      delete[] M_inv[i];
    }

  delete[] M_inv;
  delete[] res;
  delete[] p;
  delete[] denv;
  delete[] zres;
  delete[] values;
  delete[] col_in;
  delete[] row_p;
}
