#include "hermes1d.h"

// ********************************************************************

// This example solves the Poisson equation -u'' - f = 0 in
// an interval (A, B), equipped with a Dirichlet boundary
// condition on the left and a Neumann BC on the right. 

// General input:
static int N_eq = 1;
int N_elem = 30;                       // number of elements
double A = 0, B = 2*M_PI;              // domain end points
int P_init = 1;                        // initial polynomal degree

// Boundary conditions
double Val_dir_left = 0;               // Dirichlet condition left
double Val_neum_right = 0;             // Neumann condition right
                                       // (derivative at end point)

// Matrix solver
const int MATRIX_SOLVER = 1;            // 0... default (LU decomposition)
                                        // 1... UMFPACK
                                        // 2... CG (no preconditioning)
                                        // Only relevant for iterative matrix solvers:
const double MATRIX_SOLVER_TOL = 1e-7;  // Tolerance for residual in L2 norm
const int MATRIX_SOLVER_MAXITER = 150;  // Max. number of iterations

// Newton's method
double NEWTON_TOL = 1e-5;
int NEWTON_MAXITER = 150;

// Function f(x)
double f(double x) {
  return sin(x);
  //return 1;
}

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate 
  // basis functions
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left);
  printf("N_dof = %d\n", mesh->assign_dofs());

  // Register weak forms
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian_vol);
  dp->add_vector_form(0, residual_vol);
  dp->add_vector_form_surf(0, residual_surf_right, BOUNDARY_RIGHT);

  // Newton's loop
  newton(dp, mesh, MATRIX_SOLVER, MATRIX_SOLVER_TOL, MATRIX_SOLVER_MAXITER,
         NEWTON_TOL, NEWTON_MAXITER);

  // Plot the solution
  Linearizer l(mesh);
  l.plot_solution("solution.gp");

  printf("Done.\n");
  return 1;
}
