#ifndef FEM1DFUNCDEF
#define FEM1DFUNCDEF

class Neumann
{
  int elem;
  int ndnum;
  double val;
public:
  void setNuem(int, int, double);
  void applyNuem(int**, double*);
};

class Dirichlet
{
  int ndnum;
  double val;
public:
  void setDir(int, double);
  void applyDir(double**, double*, int);
};
  
void sfcn(double s, double *N, double *dNds);
void Qvf(int s, int *elem, double *node, double **B, double& Jac);
void CRS (double **A, double *values, int *col_in, int *row_p, int n, int *cs);
void CG_Solve(double **A, double *b, double *u, int n);

#endif
