1D Finite Element code which

**********************************************
Solves

  d^2u/dx^2 + u = 0
  
  u(0) = 0         u'(3.1416) = -1
**********************************************
using 1D linear lagrange shape functions and two point gaussian integration

Analytical Solution : u = sin(x)

Instructions to run:

1) Download all files in a single directory on any linux machine
2) Open the terminal in downloaded directory and type : make fem
        (An executable file named fem would be generated)
3) Type ./fem to run the program
4) Enter number of elements for finite elements
5) Results will be written in U.txt
