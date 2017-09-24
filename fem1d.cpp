#include <iostream>
#include <fstream>
#include <cassert>
#include "fem1dfunc.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  // Initialization
  double L = 3.1415926535897;
  int num_elem;
  cout << "No. of Elements = ";
  cin >> num_elem;
  
  
  // Nodal coordinates
  double *node;
  node = new double [num_elem+1];  
  for (int i=0; i<=num_elem; i++)
    {
      node[i] = i*L/num_elem;      
    }

  int num_nodes = num_elem + 1;

  // Connectivity
  int **elem;
  elem = new int* [num_elem];
  for (int i=0; i<num_elem; i++)
    {
      elem[i] = new int [2];
      elem[i][0] = i;
      elem[i][1] = i + 1;     
    }

  // Neumann BC setting
  Neumann NBC;
  NBC.setNuem(num_elem - 1, 1, -1);

  // Dirichlet BC setting
  Dirichlet DBC;
  DBC.setDir(0,0);
  
  // Gauss quadrature points and weightings
  double s[2] = {-0.577350269189626,0.577350269189626};
  int w[2] = {1,1};

   // K U F Initialization
  double **K, *F, *U;
  K = new double *[num_nodes];
  F = new double [num_nodes];
  U = new double [num_nodes];

  for (int i=0; i<num_nodes; i++)
    {
      K[i] = new double [num_nodes];
      F[i] = 0;
      U[i] = 0;
    }
  
  for (int i=0; i<num_nodes; i++)
    {
      for (int j=0; j<num_nodes; j++)
	{
	  K[i][j] = 0;
	}
    }

 
  // Stiffness Matrix Construction
  double Kl[4][4], **B, Jac;
  B = new double* [4];
  for (int i = 0; i < 4; i++)
    B[i] = new double [4];

  int jG, kG;
  
  for (int i = 0; i < num_elem; i++)
    {
      // Initialize Local matrix to zeros
      for (int j = 0; j < 2; j++)
	for (int k = 0; k < 2; k++)
	  Kl[j][k] = 0;

      for (int j = 0; j < 2; j++)
	{
	  Qvf(s[j], elem[i], node, B, Jac);

	  for (int k = 0; k < 2; k++)
	    for (int l = 0; l < 2; l++)
	      Kl[k][l] += w[j]*B[k][l]*Jac;
	}

      // Global Stiffness Matrix
      for (int j = 0; j < 2; j++)
	{
	  for (int k = 0; k < 2; k++)
	    {
	      jG = elem[i][j];
	      kG = elem[i][k];
	      K[jG][kG] += Kl[j][k];
	    }
	}      
    }

  // Applying Neumann boundary condition
  NBC.applyNuem(elem, F);

  // Dirichlet BC and making symmetric K
  DBC.applyDir(K, F, num_nodes);
  
  // Solve KU = F
  CG_Solve(K, F, U, num_nodes);

  // Writing U to a file
  ofstream write_U("U.txt");
  assert(write_U.is_open());
  for (int i = 0; i < num_nodes; i++)
    {
      write_U << node[i] << "," << U[i] << "\n";
    }
  write_U.close();
 
  // Deallocating memory
  delete[] node;

  for (int i = 0; i<num_elem; i++)
    {
      delete[] elem[i];
    }
  delete[] elem;

  for (int i=0; i<num_nodes; i++)
    {
      delete[] K[i];
    }

  delete[] K;
  delete[] U;
  delete[] F;

  for (int i = 0; i < 4; i++)
    delete[] B[i];

  delete[] B;
   
  return 0;
}
