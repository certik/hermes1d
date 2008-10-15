#ifndef cassembly_h
#define cassembly_h

#include "matrix.h"

class System
{
public:
    System()
    {
        this->A = NULL;
        this->B = NULL;
    }
    ~System()
    {
        if (this->A) delete this->A;
        if (this->B) delete this->B;
    }
    void set_mesh(double *mesh, int nmesh)
    {
        this->mesh = mesh;
        this->nmesh = nmesh;
        this->A = new SparseMatrix(nmesh*3+10);
        this->B = new SparseMatrix(nmesh*3+10);
    }
    void print_info();
    void assemble();

    double h(int i)
    {
        if (i + 1 > this->nmesh)
            return -1;
        return this->mesh[i+1] - this->mesh[i];
    }
    void set_dof_A(int i, int j, double value);
    void set_dof_B(int i, int j, double value);

    double bilinear_form_A(int i, int j);
    double bilinear_form_B(int i, int j);

    double int_grad_u_grad_v(int i, int j);
    double int_u_v(int i, int j);
    double int_grad_u_v_over_x(int i, int j);

    double *mesh;
    int nmesh;
    SparseMatrix *A, *B;
};

#endif
