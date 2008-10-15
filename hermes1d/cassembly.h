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
        return this->mesh[i+1] - this->mesh[i];
    }
    void set_dof_A(int i, int j, double value);
    void set_dof_B(int i, int j, double value);

    double *mesh;
    int nmesh;
    SparseMatrix *A, *B;
};

#endif
