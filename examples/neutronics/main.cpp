//#define NDEBUG

#include "hermes1d.h"

#include <map>
#include <cassert>

using std::cout;
using std::endl;

// ********************************************************************

// This example solves the eigenvalue problem for the neutron diffusion equation 
//	-(D.u')' + Sa.u = 1/k.nSf.u 
// in an environment composed of three slabs - inner core, outer core and a
// reflector. Reflective condition is prescribed on the left (homogeneous 
// Neumann BC) and there is vacuum on the outside of the reflector on the right
// (modelled by a Newton BC "albedo.u + D.u' = 0").

bool flag = false;				// flag for debugging purposes
bool verbose = false;

// General input:
int N_subdiv_inner = 2;      // Equidistant subdivision of the inner core macroelement
int N_subdiv_outer = 2;      // Equidistant subdivision of the outer core macroelement
int N_subdiv_reflector = 1;     // Equidistant subdivision of the reflector macroelement
int P_init_inner = 3;        // Initial polynomal degree in inner core (material 0)
int P_init_outer = 3;        // Initial polynomal degree in outer core (material 1)
int P_init_reflector = 3;       // Initial polynomal degree in reflector (material 2)
int Max_SI = 1000;           // Max. number of eigenvalue iterations
int N_SLN = 2;               // Number of solutions
double K_EFF = 1.0;          // Initial approximation

// Geometry and materials
const int N_MAT = 3;			           // Number of macroelements with different materials
const int N_GRP = 1;			           // Number of energy groups in multigroup approximation
double interfaces[N_MAT+1] = { 0, 50, 100, 125 };  // Coordinates of material regions interfaces [cm]
int Marker_inner = 0;                              // Material marker for inner core elements
int Marker_outer = 1;                              // Material marker for outer core elements
int Marker_reflector = 2;                             // Material marker for reflector elements

// Matrix solver
const int MATRIX_SOLVER = 1;            // 0... default (LU decomposition)
                                        // 1... UMFPACK
                                        // 2... CG (no preconditioning)
                                        // Only relevant for iterative matrix solvers:
const double MATRIX_SOLVER_TOL = 1e-7;  // Tolerance for residual in L2 norm
const int MATRIX_SOLVER_MAXITER = 150;  // Max. number of iterations

// Newton's method
double NEWTON_TOL = 1e-5;               // tolerance for the Newton's method
int NEWTON_MAXITER = 150;               // max. number of Newton iterations
double TOL_SI = 1e-8;                   // tol. for the source (eigenvalue) iteration

// Boundary conditions
double Val_neumann_left = 0.0;		// total reflection on the left (zero Neumann)
double Val_albedo_right = 0.5; 		// vacuum on the right

// Physical properties of each material type
static double D[N_GRP][N_MAT] = 	
	{ { 0.650, 0.750, 1.150 } };		// diffusion coefficient
static double Sa[N_GRP][N_MAT] = 	
	{ { 0.120, 0.100, 0.010 } };		// absorption cross-section
static double nSf[N_GRP][N_MAT] = 
	{ { 0.185, 0.150, 0.000 } };		// fission-yield cross section (\nu \Sigma_f)
static double chi[N_GRP] = 
	{ 1.0 };				// fission spectrum (for multigroup calc.)

// Other physical properties
static double nu = 2.43; 		// Nean number of neutrons released by fission
static double eps = 3.204e-11;		// Nean energy release of each fission evt [J]

/******************************************************************************/

// Calculate \int \nu \Sigma_f(x) u(x) over an element 'e'
double calc_elem_fission_yield(Element *e)
{
  // get solution values at quadrature points
  double val_phys[MAX_EQN_NUM][MAX_QUAD_PTS_NUM];
  double der_phys[MAX_EQN_NUM][MAX_QUAD_PTS_NUM];
  int order = e->p;   // this is enough since "nSf" is constant in elements
  e->get_solution_quad(0,  order, val_phys, der_phys);
  // get quadrature weights
  double phys_x[MAX_QUAD_PTS_NUM];
  double phys_weights[MAX_QUAD_PTS_NUM];
  int pts_num;
  create_phys_element_quadrature(e->x1, e->x2, order, phys_x, phys_weights, &pts_num); 
  // numerical quadrature in element 'e'
  int n_grp = e->n_eq;
  double yield = 0;
  for (int i = 0; i < pts_num; i++) {
    double val = 0;
    int m = e->marker;
    //if (flag) cout << "\t" << phys_x[i] << "\t" << m << endl; 
    for (int g = 0; g < n_grp; g++) val += nSf[g][m] * val_phys[g][i];
    yield += val * phys_weights[i];
  }
  return yield;
}

// Calculate \int_\Omega \nu \Sigma_f(x) u(x) over the entire mesh
double calc_fission_yield(Mesh* mesh)
{
  double fis_yield = 0;
  Iterator *I = new Iterator(mesh);
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
    //if (flag) printf("%d : (%f,%f) \n", e->id, e->x1, e->x2);
    //if (flag) cout << e->id << " : (" << e->x1 << "," << e->x2 << ")" << endl;
    //if (flag) cout << "hi";
    fis_yield += calc_elem_fission_yield(e);
  }
  delete I;
  return fis_yield;
}

// Normalize the eigenfunction representing the neutron flux so that the total
// power it generates equals to 'desired_power' [W]
void normalize_to_power(Mesh* mesh, double desired_power)
{
  // Calculate total power generated by the computed flux 'u': 
  // P(u) = \eps \int_\Omega \Sigma_f(x) u(x)
  double P = eps * calc_fission_yield(mesh) / nu;
	
  // Calculate normalization constant 'c', so that P(c u) = 'desired_power' 
  double c = desired_power / P;
	
  // Multiply the computed flux by the normalization constant
  multiply_dofs_with_constant(mesh, c);
}

// Weak forms
#include "forms.cpp"

/******************************************************************************/

int main() {
  // Three macroelements are defined above via the interfaces[] array
  // poly_orders[]... initial poly degrees of macroelements
  // material_markers[]... material markers of macroelements
  // subdivisions[]... equidistant subdivision of macroelements
  int poly_orders[N_MAT] = {P_init_inner, P_init_outer, P_init_reflector };
  int material_markers[N_MAT] = {Marker_inner, Marker_outer, Marker_reflector };
  int subdivisions[N_MAT] = {N_subdiv_inner, N_subdiv_outer, N_subdiv_reflector };

  // Create coarse mesh, enumerate basis functions
  Mesh *mesh = new Mesh(N_MAT, interfaces, poly_orders, material_markers, subdivisions, N_GRP, N_SLN);
  printf("N_dof = %d\n", mesh->assign_dofs());

  // Initial approximation: u = 1
  double K_EFF_old;
  double init_val = 1.0;
  set_vertex_dofs_constant(mesh, init_val, 0);
  
  // Register weak forms
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian_vol_inner, Marker_inner);
  dp->add_matrix_form(0, 0, jacobian_vol_outer, Marker_outer);
  dp->add_matrix_form(0, 0, jacobian_vol_reflector, Marker_reflector);
  dp->add_vector_form(0, residual_vol_inner, Marker_inner);
  dp->add_vector_form(0, residual_vol_outer, Marker_outer);
  dp->add_vector_form(0, residual_vol_reflector, Marker_reflector);
  dp->add_vector_form_surf(0, residual_surf_left, BOUNDARY_LEFT);
  dp->add_matrix_form_surf(0, 0, jacobian_surf_right, BOUNDARY_RIGHT);
  dp->add_vector_form_surf(0, residual_surf_right, BOUNDARY_RIGHT);

  // Source iteration (power method)
  for (int i = 0; i < Max_SI; i++)
  {	
    // Obtain fission source
    int current_solution = 0, previous_solution = 1;
    copy_dofs(current_solution, previous_solution, mesh);
		
    // Newton's loop		
    newton(dp, mesh, MATRIX_SOLVER, MATRIX_SOLVER_TOL, MATRIX_SOLVER_MAXITER,
           NEWTON_TOL, NEWTON_MAXITER, verbose);
	
    // Update the eigenvalue
    K_EFF_old = K_EFF;
    K_EFF = calc_fission_yield(mesh);		
    if (verbose) printf("K_EFF_%d = %f\n", i, K_EFF);
		
    if (fabs(K_EFF - K_EFF_old)/K_EFF < TOL_SI) break;
  }
	
  flag = true;
  
  // Plot the critical (i.e. steady-state) neutron flux
  Linearizer l(mesh);
  l.plot_solution("solution.gp");
  
  // Normalize so that the absolute neutron flux generates 320 Watts of energy
  // (note that, using the symmetry condition at the origin, we've solved for  
  // flux only in the right half of the reactor) 
  normalize_to_power(mesh, 320/2.);	

  // Plot the solution and mesh
  l.plot_solution("solution_320W.gp");	
  mesh->plot("mesh.gp");

  printf("K_EFF = %f\n", K_EFF);

  printf("Done.\n");
  return 1;
}
