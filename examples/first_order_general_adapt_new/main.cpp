#include "hermes1d.h"

// ********************************************************************

// This example uses automatic hp-adaptivity to solve the general 
// first-order equation y' = f(y, x) in an interval (A, B), equipped 
// with the initial condition y(A) = YA. The function f can be linear
// or nonlinear in 'y', as long as it is differentiable
// with respect to this variable (needed for the Newton's method). 

// General input:
const int N_eq = 1;                     // Number of equations
const int N_elem = 5;                   // Number of elements
const double A = 0, B = 10;             // Domain end points
const double YA = 1;                    // Equation parameter
const int P_init = 1;                   // Initial polynomal degree

// Stopping criteria for Newton
const double TOL_NEWTON_COARSE = 1e-8;  // Coarse mesh
const double TOL_NEWTON_REF = 1e-8;     // Fine mesh

// Adaptivity
const int ADAPT_TYPE = 0;         // 0... hp-adaptivity
                                  // 1... h-adaptivity
                                  // 2... p-adaptivity
const double THRESHOLD = 0.7;     // Refined will be all elements whose error
                                  // is greater than THRESHOLD*max_elem_error
const double TOL_ERR_REL = 1e-4;  // Tolerance for the relative error between 
                                  // the coarse and fine mesh solutions
const int NORM = 0;               // To measure errors:
                                  // 1... H1 norm
                                  // 0... L2 norm
 
// Right-hand side function f(y, x)
double f(double y, double x) {
  //return -y; // with y(0)=1, exact solution is y=exp(-x)
  return -y*y; // with y(0)=1, exact solution is y=1/(x+1)
}

// y-derivative of dfdy(y, x)
double dfdy(double y, double x) {
  //return -1;
  return -2*y;
}

// Exact solution
const int EXACT_SOL_PROVIDED = 1;
double exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) {
  u[0] = 1./(x+1);
  dudx[0] = -1/((x+1)*(x+1));
}

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate 
  // basis functions
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, YA);
  printf("N_dof = %d\n", mesh->assign_dofs());

  // Create discrete problem on coarse mesh
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian);
  dp->add_vector_form(0, residual);

  // Convergence graph wrt. the number of degrees of freedom
  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Convergence History", "Degrees of Freedom", "Error [%]");
  graph.add_row("exact error", "k", "-", "o");
  graph.add_row("error estimate", "k", "--");

  // Main adaptivity loop
  int adapt_iterations = 1;
  Mesh *mesh_ref;
  while(1) {
    printf("============ Adaptivity step %d ============\n", adapt_iterations); 

    // Newton's loop on coarse mesh
    int success, iter_num;
    success = newton(dp, mesh, TOL_NEWTON_COARSE, iter_num);
    if (!success) error("Newton's method did not converge."); 
    printf("Finished initial coarse mesh Newton's iteration (%d iter).\n", 
           iter_num);

    /***** Identify elements whose refinement causes biggest change in
           the solution (later to be replaced with quantity of interest )*****/

    // For every element perform fast trial refinement
    // and calculate an array of the corresponding errors
    int n_elem = mesh->get_n_active_elem();
    double ftr_errors[MAX_ELEM_NUM];
    for (int i=0; i < n_elem; i++) {
      mesh_ref = mesh->replicate();
      // refine one element starting with element 'i'
      mesh_ref->reference_refinement(i, 1); 
      printf("Elem [%d]: fine mesh created (%d DOF).\n", i, 
                 mesh_ref->assign_dofs());

      // Transfer solution from 'mesh' to 'mesh_ref' 
      // FIXME: this should be done on one element only,
      //        in all others the solution remains the same !!! 
      transfer_solution_forward(mesh, mesh_ref);
      printf("Elem [%d]: coarse mesh solution copied to fine mesh.\n", i);

      // Newton's loop on fine mesh
      success = newton(dp, mesh_ref, TOL_NEWTON_REF, iter_num);
      if (!success) error("Newton's method did not converge."); 
      printf("Elem [%d]: finished fine mesh Newton's iteration (%d iter).\n", 
             i, iter_num);

      // Print fine solution (enumerated) 
      Linearizer *lxx = new Linearizer(mesh_ref);
      char out_filename[255];
      sprintf(out_filename, "solution_ref_%d.gp", i);
      lxx->plot_solution(out_filename);
      delete lxx;

      // Calculate norm of the difference between the locally refined 
      // and coarse mesh solutions.
      double err_est_array[MAX_ELEM_NUM];
      ftr_errors[i] = calc_elem_est_errors(NORM, mesh, mesh_ref, 
                                           err_est_array);

      delete mesh_ref;
    }

    /***** Select the best refinement for each element to be refined *****/

    // Create globally refined mesh (all elements in 'h' and 'p')
    mesh_ref = mesh->replicate();
    int start_elem_id = 0; 
    int num_to_ref = mesh_ref->get_n_active_elem();
    mesh_ref->reference_refinement(start_elem_id, num_to_ref);
    printf("Fine mesh created (%d DOF).\n", mesh_ref->get_n_dof());

    // Transfer coarse mesh solution to the fine mesh
    transfer_solution_forward(mesh, mesh_ref);
    printf("Coarse mesh solution copied to fine mesh.\n");

    // Newton's loop on fine mesh
    success = newton(dp, mesh_ref, TOL_NEWTON_REF, iter_num);
    if (!success) error("Newton's method did not converge."); 
    printf("Finished fine mesh Newton's iteration (%d iter).\n", 
           iter_num);

    // Estimate element errors based on the difference between the 
    // fine mesh and coarse mesh solutions. 
    double err_est_array[MAX_ELEM_NUM]; 
    double err_est_total = calc_elem_est_errors(NORM, mesh, mesh_ref, 
                                                err_est_array);

    // Calculate the norm of the fine mesh solution
    double ref_sol_norm = calc_approx_sol_norm(NORM, mesh_ref);

    // Calculate an estimate of the global relative error
    double err_est_rel = err_est_total/ref_sol_norm;
    printf("Relative error (est) = %g %%\n", 100.*err_est_rel);

    // If exact solution available, also calculate exact error
    if (EXACT_SOL_PROVIDED) {
      // Calculate element errors wrt. exact solution
      double err_exact_total = calc_exact_sol_error(NORM, mesh, exact_sol);
     
      // Calculate the norm of the exact solution
      // (using a fine subdivision and high-order quadrature)
      int subdivision = 500; // heuristic parameter
      int order = 20;        // heuristic parameter
      double exact_sol_norm = calc_exact_sol_norm(NORM, exact_sol, N_eq, A, B,
                                                  subdivision, order);
      // Calculate an estimate of the global relative error
      double err_exact_rel = err_exact_total/exact_sol_norm;
      printf("Relative error (exact) = %g %%\n", 100.*err_exact_rel);
      graph.add_values(0, mesh->get_n_dof(), 100 * err_exact_rel);
    }

    // add entry to DOF convergence graph
    graph.add_values(1, mesh->get_n_dof(), 100 * err_est_rel);

     // Decide whether the relative error is sufficiently small
    if(err_est_rel*100 < TOL_ERR_REL) break;

    // debug
    //if (adapt_iterations == 3) break;

    /***** Perform the refinements *****/

    // Returns updated coarse and fine meshes, with the last 
    // coarse and fine mesh solutions on them, respectively. 
    // The coefficient vectors and numbers of degrees of freedom 
    // on both meshes are also updated. 
    adapt(NORM, ADAPT_TYPE, THRESHOLD, ftr_errors,
          mesh, mesh_ref);

    adapt_iterations++;
  }

  // Plot meshes, results, and errors
  adapt_plotting(mesh, mesh_ref,
                 NORM, EXACT_SOL_PROVIDED, exact_sol);

  // Save convergence graph
  graph.save("conv_dof.gp");

  printf("Done.\n");
  return 1;
}

















