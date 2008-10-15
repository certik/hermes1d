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
        this->A = new SparseMatrix(nmesh-1);
        this->B = new SparseMatrix(nmesh-1);
    }
    void print_info();
    void assemble();

    double h(int i)
    {
        return this->mesh[i+1] - this->mesh[i];
    }

    double *mesh;
    int nmesh;
    SparseMatrix *A, *B;
};

#endif
